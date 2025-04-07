package main

import (
	"flag"
	"fmt"
	"genomics/comparison"
	"genomics/database"
	"genomics/genomes"
	"genomics/pileup"
	"genomics/utils"
	"log"
	"os"
	"path"
	"slices"
	"strings"
	"time"
)

/*
Find how many locations have two or more alleles with a depth of more than
minDepth
*/
func FindQS(pu *pileup.Pileup, minDepth int) []utils.OneBasedPos {
	ret := make([]utils.OneBasedPos, 0)
outer:
	for i := 0; i <= pu.MaxPos; i++ {
		record := pu.Get(i)
		if record == nil {
			continue
		}
		count := 0
		for _, read := range record.Reads {
			if read.Depth >= minDepth {
				count++
			}
			if count >= 2 {
				ret = append(ret, utils.OneBasedPos(record.Pos+1))
				continue outer
			}
		}
	}
	return ret
}

func FindPileup(record *database.Record, root string) *pileup.Pileup {
	path := path.Join(root, fmt.Sprintf("%s-WH1-index.txt.gz", record.SRAs()))
	if _, err := os.Stat(path); err != nil {
		// log.Printf("%s doesn't exist\n", path)
		return nil
	}
	ret, err := pileup.Parse2(path)
	if err != nil {
		log.Fatal(err)
	}
	return ret
}

func LoadRecords(db *database.Database, fname string) []database.Id {
	ret := make([]database.Id, 0)
	utils.Lines(fname, func(line string, lineErr error) bool {
		ids := db.GetByAccession(line)
		ret = append(ret, ids...)
		return true
	})
	return ret
}

type ReadHandler interface {
	Process(*database.Record, *pileup.Pileup)
}

type Allele struct {
	pos int
	nt  byte
}

type CountAll struct {
	ref      *genomes.Genomes
	minDepth int
	counts   map[Allele]int
	dates    map[Allele][]time.Time
	cutoff   time.Time
}

func (c *CountAll) Init(ref *genomes.Genomes, minDepth int, cutoff time.Time) {
	c.minDepth = minDepth
	c.counts = make(map[Allele]int)
	c.dates = make(map[Allele][]time.Time)
	c.ref = ref
}

func (c *CountAll) Process(record *database.Record, pu *pileup.Pileup) {
	// We're looking for anywhere we have reads exceeding our min-depth that
	// differ from WH1.
	for pos := 0; pos <= pu.MaxPos; pos++ {
		pur := pu.Get(pos)
		if pur == nil {
			continue
		}
		for _, read := range pur.Reads {
			if read.Depth < c.minDepth {
				continue
			}
			if c.ref.Nts[0][pos] == read.Nt {
				continue
			}
			allele := Allele{pos, read.Nt}
			c.counts[allele]++

			dates, there := c.dates[allele]
			if !there {
				dates = make([]time.Time, 0, 1)
			}
			dates = append(dates, record.CollectionDate)
			c.dates[allele] = dates
		}
	}
}

func (c *CountAll) SortDates() {
	for _, v := range c.dates {
		slices.SortFunc(v, func(a, b time.Time) int {
			return a.Compare(b)
		})
	}
}

func dateString(dates []time.Time) string {
	return "-" // comment this out if you really want the dates
	s := make([]string, len(dates))
	for i, date := range dates {
		s[i] = fmt.Sprintf("%d",
			int(date.Sub(utils.Date(2020, 1, 1)).Hours()/24))
	}
	return strings.Join(s, " ")
}

func (c *CountAll) Display() {
	og := NewOutgroup()

	type result struct {
		Allele
		silent bool
		alts   Alts
		count  int
		dates  []time.Time
	}
	results := make([]result, 0, len(c.counts))

	for allele, count := range c.counts {
		silent, _, _ := genomes.IsSilentWithReplacement(c.ref,
			allele.pos, 0, 0, []byte{allele.nt})
		alts := og.Get(allele.pos)
		results = append(results,
			result{allele, silent, alts, count, c.dates[allele]})
	}

	slices.SortFunc(results, func(a, b result) int {
		if a.count > b.count {
			return -1
		}
		if a.count < b.count {
			return 1
		}
		return 0
	})

	for _, r := range results {
		var aaChange string
		var silent string
		refNt := c.ref.Nts[0][r.pos]

		if !r.silent {
			_, oldProt, newProt, err := genomes.ProteinChange(c.ref, r.pos,
				0, 0, []byte{r.nt})
			if err == nil {
				mut := comparison.Mut{oldProt[0], newProt[0], r.pos}
				aaChange = mut.ToString(c.ref.Orfs)
			}
		} else {
			silent = "*"
		}
		fmt.Printf("%c%d%c%s %s %s %d on %s\n", refNt, r.pos+1,
			r.nt, silent, aaChange,
			r.alts.ToString(), r.count, dateString(r.dates))
	}
}

type Display struct {
	minDepth int
}

func (d *Display) Init(minDepth int) {
	d.minDepth = minDepth
}

func (d *Display) Process(record *database.Record, pu *pileup.Pileup) {
	contents := Classify(pu, d.minDepth)
	qsLocations := FindQS(pu, contents.QSDepth())

	fmt.Println(record.CollectionDate.Format(time.DateOnly),
		record.GisaidAccession, record.SRAs(),
		record.Region, record.Country,
		contents.ToString(), len(qsLocations))
}

func ProcessReads(db *database.Database,
	ids []database.Id, prefix string,
	handler ReadHandler) {

	root := "/fs/bowser/genomes/raw_reads/"
	for _, id := range ids {
		record := db.Get(id)
		pu := FindPileup(record, path.Join(root, prefix))
		if pu == nil {
			continue
		}
		handler.Process(record, pu)
	}
}

func DisplayReads(db *database.Database,
	ids []database.Id, minDepth int, prefix string) {
	fmt.Println("Date AccNo SRA Region Country class 8782:C|T 28144:C|T QSlocs")

	var d Display
	d.Init(minDepth)
	ProcessReads(db, ids, prefix, &d)
}

func CountEverything(db *database.Database,
	ids []database.Id, minDepth int, prefix string, cutoff time.Time) {
	var c CountAll
	ref := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	c.Init(ref, minDepth, cutoff)
	ProcessReads(db, ids, prefix, &c)
	c.SortDates()
	c.Display()
}

/*
Find the apparent CC and TT sequences in the database (GISAID 2020) which have
reads available.
*/
func FindSequences(db *database.Database, showClass string) {
	counts := make(map[string]int)
	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}

		if len(r.SRA) == 0 {
			return false
		}

		C1, C2 := true, false
		for _, m := range r.NucleotideChanges {
			if m.Pos == 8782 && m.To == 'T' {
				C1 = false
			}
			if m.Pos == 28144 && m.To == 'C' {
				C2 = true
			}
		}

		var class string

		display := func() {
			fmt.Println(class, r.SRAs(), len(r.NucleotideChanges),
				r.GisaidAccession, r.Country,
				r.CollectionDate.Format(time.DateOnly))
		}

		if C1 {
			if C2 {
				class = "CC"
			} else {
				class = "CT" // aka Lin B
			}
		} else {
			if C2 {
				class = "TC" // aka Lin A
			} else {
				class = "TT"
			}
		}

		display()
		/*
			if class == showClass {
				display()
			}
		*/

		counts[class]++
		return false
	})

	fmt.Println(counts)
}

func ShowLocations(db *database.Database, ids []database.Id) {
	countries := make(map[string]int)
	db.Filter(utils.ToSet(ids), func(r *database.Record) bool {
		countries[r.Country]++
		return true
	})

	type result struct {
		country string
		count   int
	}
	results := make([]result, 0, len(countries))
	for k, v := range countries {
		results = append(results, result{k, v})
	}
	slices.SortFunc(results, func(a, b result) int {
		if a.count > b.count {
			return -1
		}
		if a.count < b.count {
			return 1
		}
		return 0
	})

	for _, v := range results {
		fmt.Println(v.country, v.count)
	}
}

func main() {
	var (
		minDepth      int
		findSequences bool
		class         string
		count         bool
		countCT       bool
		cutoffS       string
		cutoff        time.Time
		locations     bool
		showPoss      bool
	)

	flag.BoolVar(&findSequences, "f", false, "Find the sequences")
	flag.StringVar(&class, "class", "first100", "Classes to show")
	flag.IntVar(&minDepth, "min-depth", 3, "Min depth")
	flag.BoolVar(&count, "countCT", false, "Count silent CT")
	flag.BoolVar(&count, "count", false, "Count everything")
	flag.StringVar(&cutoffS, "cutoff", "", "Count before date")
	flag.BoolVar(&locations, "locations", false, "Just show locations")
	flag.BoolVar(&showPoss, "show-poss", false, "Show possible silent")
	flag.Parse()

	if cutoffS != "" {
		var err error
		cutoff, err = time.Parse(time.DateOnly, cutoffS)
		if err != nil {
			log.Fatal(err)
		}
	}

	if showPoss {
		ShowPossible()
		return
	}

	db := database.NewDatabase()

	if findSequences {
		FindSequences(db, class)
		return
	}

	records := LoadRecords(db, class)
	if locations {
		ShowLocations(db, records)
	} else if countCT {
		CountReadsCT(db, records, minDepth, class, cutoff)
	} else if count {
		CountEverything(db, records, minDepth, class, cutoff)
	} else {
		DisplayReads(db, records, minDepth, class)
	}
}
