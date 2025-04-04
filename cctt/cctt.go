package main

import (
	"flag"
	"fmt"
	"genomics/comparison"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/pileup"
	"genomics/stats"
	"genomics/utils"
	"log"
	"os"
	"path"
	"slices"
	"time"
)

type Contents struct {
	C8782 int // The depth of C at 8782
	T8782 int

	C28144 int
	T28144 int

	Classification string
	Significance   float64
}

func (c *Contents) ToString() string {
	return fmt.Sprintf("%s %d|%d %d|%d %.4g", c.Classification,
		c.C8782, c.T8782, c.C28144, c.T28144, c.Significance)
}

/*
Suppose you find much more C than T in both locations. What is the probability
of doing that under the assumption that what you have is a mix of CT and TC? If
it were the latter, you would expect C/T to roughly equal T/C in the two
locations.
*/
func (c *Contents) CalcSignificance() {
	var ct stats.ContingencyTable
	ct.Init(c.C8782, c.T8782, c.T28144, c.C28144)
	ct.FisherExact(stats.TWO_SIDED)
	c.Significance = ct.P
}

func Classify(pu *pileup.Pileup, minDepth int) Contents {
	pos8782 := pu.Get(8782 - 1)
	pos28144 := pu.Get(28144 - 1)

	ret := Contents{
		pos8782.GetDepthOf('C'), pos8782.GetDepthOf('T'),
		pos28144.GetDepthOf('C'), pos28144.GetDepthOf('T'),
		"-", 1.0,
	}

	if ret.C8782 >= minDepth &&
		ret.T8782 < minDepth &&
		ret.C28144 >= minDepth &&
		ret.T28144 < minDepth {
		ret.Classification = "CC*"
	} else if ret.T8782 >= minDepth &&
		ret.C8782 < minDepth &&
		ret.T28144 >= minDepth &&
		ret.C28144 < minDepth {
		ret.Classification = "TT*"
	} else if ret.C8782 > ret.T8782 &&
		ret.C28144 > ret.T28144 &&
		ret.C8782 >= minDepth &&
		ret.C28144 >= minDepth {
		ret.Classification = "CC>"
	} else if ret.T8782 > ret.C8782 &&
		ret.T28144 > ret.C28144 &&
		ret.T8782 >= minDepth &&
		ret.T28144 >= minDepth {
		ret.Classification = "TT>"
	}
	ret.CalcSignificance()
	return ret
}

func (c *Contents) QSDepth() int {
	return utils.Max(utils.Min(c.C8782, c.T8782), utils.Min(c.C28144, c.T28144))
}

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

type CountAll struct {
	minDepth int
	counts   map[int]int // How many in each position
	silentCT []int       // All possible silent C->T positions
	cutoff   time.Time
	byDate   map[int]int
}

type Allele struct {
	pos int
	nt  byte
}

type CountAll2 struct {
	ref      *genomes.Genomes
	minDepth int
	counts   map[Allele]int
	cutoff   time.Time
}

func (c *CountAll) Init(ref *genomes.Genomes, minDepth int, cutoff time.Time) {
	c.minDepth = minDepth
	c.cutoff = cutoff
	muts := mutations.PossibleSilentMuts(ref, 0)
	c.silentCT = make([]int, 0)
	for _, mut := range muts {
		if mut.From == 'C' && mut.To == 'T' {
			c.silentCT = append(c.silentCT, mut.Pos)
		}
	}
	c.counts = make(map[int]int)
	c.byDate = make(map[int]int)
}

func (c *CountAll2) Init(ref *genomes.Genomes, minDepth int, cutoff time.Time) {
	c.minDepth = minDepth
	c.counts = make(map[Allele]int)
	c.ref = ref
}

func (c *CountAll) Process(record *database.Record, pu *pileup.Pileup) {
	for _, pos := range c.silentCT {
		pur := pu.Get(pos)
		if pur == nil {
			continue
		}
		if !c.cutoff.IsZero() {
			if record.CollectionDate.Compare(c.cutoff) > 0 {
				continue
			}
		}
		for _, read := range pur.Reads {
			if read.Depth >= c.minDepth && read.Nt == 'T' {
				c.counts[pur.Pos] += read.Depth
				delta := int(record.CollectionDate.Sub(utils.Date(2020,
					1, 1)).Hours() / 24)
				c.byDate[delta] = len(c.counts)
			}
		}
	}

}

func (c *CountAll2) Process(record *database.Record, pu *pileup.Pileup) {
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
		}
	}
}

func (c *CountAll) Display() {
	og := NewOutgroup()

	type count struct {
		pos, count int
	}
	counts := make([]count, 0, len(c.counts))
	for k, v := range c.counts {
		counts = append(counts, count{k, v})
	}
	slices.SortFunc(counts, func(a, b count) int {
		if a.count > b.count {
			return -1
		}
		if a.count < b.count {
			return 1
		}
		return 0
	})
	fmt.Printf("%d possible silent C->Ts\n", len(c.silentCT))
	for _, count := range counts {
		fmt.Printf("%d %d %s\n", count.pos+1,
			count.count, og.Get(count.pos).ToString())
	}

	fmt.Println("By Date")

	type dateCount struct {
		days int // since 2020-01-01
		seen int // number of different C->T muts seen
	}

	dateCounts := make([]dateCount, 0, len(c.byDate))
	for k, v := range c.byDate {
		dateCounts = append(dateCounts, dateCount{k, v})
	}
	slices.SortFunc(dateCounts, func(a, b dateCount) int {
		if a.days < b.days {
			return -1
		}
		if a.days > b.days {
			return 1
		}
		return 0
	})
	for _, dc := range dateCounts {
		fmt.Println(dc.days, dc.seen)
	}
}

func (c *CountAll2) Display() {
	og := NewOutgroup()

	type result struct {
		Allele
		silent bool
		alts   Alts
		count  int
	}
	results := make([]result, 0, len(c.counts))

	for allele, count := range c.counts {
		silent, _, _ := genomes.IsSilentWithReplacement(c.ref,
			allele.pos, 0, 0, []byte{allele.nt})
		alts := og.Get(allele.pos)
		results = append(results, result{allele, silent, alts, count})
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
		fromNt := c.ref.Nts[0][r.pos]

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
		fmt.Printf("%c%d%c%s %s %s %d\n", fromNt, r.pos+1,
			r.nt, silent, aaChange, r.alts.ToString(), r.count)
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

func CountReadsCT(db *database.Database,
	ids []database.Id, minDepth int, prefix string, cutoff time.Time) {
	var c CountAll
	ref := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	c.Init(ref, minDepth, cutoff)
	ProcessReads(db, ids, prefix, &c)
	c.Display()
}

func CountEverything(db *database.Database,
	ids []database.Id, minDepth int, prefix string, cutoff time.Time) {
	var c CountAll2
	ref := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	c.Init(ref, minDepth, cutoff)
	ProcessReads(db, ids, prefix, &c)
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
	)

	flag.BoolVar(&findSequences, "f", false, "Find the sequences")
	flag.StringVar(&class, "class", "TT", "Classes to show")
	flag.IntVar(&minDepth, "min-depth", 3, "Min depth")
	flag.BoolVar(&count, "countCT", false, "Count silent CT")
	flag.BoolVar(&count, "count", false, "Count everything")
	flag.StringVar(&cutoffS, "cutoff", "", "Count before date")
	flag.BoolVar(&locations, "locations", false, "Just show locations")
	flag.Parse()

	if cutoffS != "" {
		var err error
		cutoff, err = time.Parse(time.DateOnly, cutoffS)
		if err != nil {
			log.Fatal(err)
		}
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
