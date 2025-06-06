package main

import (
	"flag"
	"fmt"
	"genomics/comparison"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/pileup"
	"genomics/utils"
	"log"
	"os"
	"path"
	"path/filepath"
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

type Count struct {
	count         int
	maxDepth      int
	maxDepthRatio float64
}

type CountAll struct {
	ref      *genomes.Genomes
	minDepth int
	counts   map[Allele]Count
	dates    map[Allele][]time.Time
	cutoff   time.Time
	nonMaj   bool

	// How often each possible silent mutation actually appears above minDepth
	possibleSilent map[mutations.Mutation]Count
	possibleNS     map[mutations.Mutation]Count

	// For each day since the start, the mean depth ratio.
	evolution map[int]Mean
}

func (c *CountAll) Init(ref *genomes.Genomes,
	minDepth int, cutoff time.Time, nonMaj bool) {
	c.minDepth = minDepth
	c.counts = make(map[Allele]Count)
	c.dates = make(map[Allele][]time.Time)
	c.evolution = make(map[int]Mean)
	c.ref = ref
	c.nonMaj = nonMaj

	c.possibleSilent = make(map[mutations.Mutation]Count)
	for _, mut := range mutations.PossibleSilentMuts(ref, 0) {
		c.possibleSilent[mut] = Count{}
	}

	c.possibleNS = make(map[mutations.Mutation]Count)
	for _, mut := range mutations.PossibleNonSilentMuts(ref, 0) {
		c.possibleNS[mut] = Count{}
	}
}

func (c *CountAll) PossibleCTTC() (int, int) {
	var ct, tc int
	for k, _ := range c.possibleSilent {
		if k.From == 'C' && k.To == 'T' {
			ct++
		} else if k.From == 'T' && k.To == 'C' {
			tc++
		}
	}
	return ct, tc
}

func (c *CountAll) Process(record *database.Record, pu *pileup.Pileup) {
	for pos := 0; pos <= pu.MaxPos; pos++ {
		pur := pu.Get(pos)
		if pur == nil {
			continue
		}
		for i, read := range pur.Reads {
			if c.nonMaj && i == 0 {
				continue
			}
			if read.Depth < c.minDepth {
				continue
			}
			if c.ref.Nts[0][pos] == read.Nt {
				continue
			}

			/*
			if pos == 21710 && read.Nt == 'T' {
				days := int(record.CollectionDate.Sub(
					utils.Date(2020, 1, 1)).Hours() / 24)
				fmt.Println("T21711", days, read.Depth)
			}
			*/

			allele := Allele{pos, read.Nt}

			depthRatio := float64(read.Depth) / float64(pur.TotalDepth)
			existing := c.counts[allele]
			c.counts[allele] = Count{existing.count + 1,
				utils.Max(existing.maxDepth, read.Depth),
				utils.Max(existing.maxDepthRatio, depthRatio)}

			dates, there := c.dates[allele]
			if !there {
				dates = make([]time.Time, 0, 1)
			}

			if record != nil {
				dates = append(dates, record.CollectionDate)
				c.dates[allele] = dates

				days := int(record.CollectionDate.Sub(
					utils.Date(2020, 1, 1)).Hours() / 24)
				mean := c.evolution[days]
				mean.Add(depthRatio)
				c.evolution[days] = mean
			}

			mut := mutations.Mutation{mutations.BaseMutation{pos, true},
				c.ref.Nts[0][pos], read.Nt}
			if _, there := c.possibleSilent[mut]; there {
				existing := c.possibleSilent[mut]
				c.possibleSilent[mut] = Count{existing.count + 1,
					utils.Max(existing.maxDepth, read.Depth),
					utils.Max(existing.maxDepthRatio, depthRatio)}
			}
			mut.Silent = false
			if _, there := c.possibleNS[mut]; there {
				existing := c.possibleNS[mut]
				c.possibleNS[mut] = Count{existing.count + 1,
					utils.Max(existing.maxDepth, read.Depth),
					utils.Max(existing.maxDepthRatio, depthRatio)}
			}
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
	return "" // comment this out if you really want the dates
	s := make([]string, len(dates))
	for i, date := range dates {
		s[i] = fmt.Sprintf("%d",
			int(date.Sub(utils.Date(2020, 1, 1)).Hours()/24))
	}
	return strings.Join(s, " ")
}

type Transition struct {
	from, to byte
}

type TransitionCount struct {
	possible      int // How many of these are possible
	total         int // total number of samples with this transition
	maxDepthRatio float64
}

func (c *CountAll) TransitionRatios() {
	transitions := make(map[Transition]TransitionCount)
	for k, v := range c.possibleSilent {
		transition := Transition{k.From, k.To}
		existing := transitions[transition]
		transitions[transition] = TransitionCount{existing.possible + 1,
			existing.total + v.count,
			utils.Max(existing.maxDepthRatio, v.maxDepthRatio)}
	}

	type result struct {
		trans Transition
		count TransitionCount
	}
	results := make([]result, 0, len(transitions))
	for k, v := range transitions {
		results = append(results, result{k, v})
	}
	utils.SortByKey(results, false, func(r result) int {
		return r.count.total
	})
	for _, r := range results {
		ratio := float64(r.count.total) / float64(r.count.possible)
		special := (r.trans.from == 'C' && r.trans.to == 'T') ||
			(r.trans.from == 'T' && r.trans.to == 'C')
		var highlight string
		if special {
			highlight = "*** "
		}
		fmt.Printf("%s%c->%c %d %d %.2f %f\n", highlight,
			r.trans.from, r.trans.to,
			r.count.possible, r.count.total,
			ratio, r.count.maxDepthRatio)
	}
}

func (c *CountAll) Display(minSamples int) {
	og := NewOutgroup()

	type result struct {
		Allele
		silent bool
		alts   Alts
		count  Count
		dates  []time.Time
	}
	results := make([]result, 0, len(c.counts))

	var revCT, revTC int // reversions to OG
	var divCT, divTC int // "diversions": mutations away from OG

	for allele, count := range c.counts {
		silent, _, _ := genomes.IsSilentWithReplacement(c.ref,
			allele.pos, 0, 0, []byte{allele.nt})
		alts := og.Get(allele.pos)
		results = append(results,
			result{allele, silent, alts, count, c.dates[allele]})
	}

	utils.SortByKey(results, false,
		func(a result) float64 { return a.count.maxDepthRatio })

	var ogRevs int
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
		ds := dateString(r.dates)
		if ds != "" {
			ds = fmt.Sprintf("on %s", ds)
		}
		var ogReversion string
		if r.nt == r.alts[0].nt {
			ogReversion = " OG"
			ogRevs++

			if r.silent {
				if refNt == 'C' && r.nt == 'T' {
					revCT++
				} else if refNt == 'T' && r.nt == 'C' {
					revTC++
				}
			}
		}

		if refNt == r.alts[0].nt {
			if r.silent {
				if refNt == 'C' && r.nt == 'T' {
					divCT++
				} else if refNt == 'T' && r.nt == 'C' {
					divTC++
				}
			}
		}

		if r.count.count >= minSamples {
			fmt.Printf("%c%d%c%s %s %s %d@%d%s (%.2f)%s\n",
				refNt, r.pos+1,
				r.nt, silent, aaChange,
				r.alts.ToString(),
				r.count.count, r.count.maxDepth, ds, r.count.maxDepthRatio,
				ogReversion)
		}
	}
	fmt.Printf("%d OG reversions out of %d (%.4f)\n", ogRevs, len(results),
		float64(ogRevs)/float64(len(results)))

	possCT, possTC := c.PossibleCTTC()
	fmt.Printf("%d/%d silent CT reversions, %d/%d silent TC reversions\n",
		revCT, possCT, revTC, possTC)
	fmt.Printf("%d/%d (%.2f) silent CT diversions, "+
		"%d/%d (%.2f) silent TC diversions\n",
		divCT, possCT,
		float64(divCT)/float64(possCT),
		divTC, possTC,
		float64(divTC)/float64(possTC))
}

func (c *CountAll) DisplayPossible(minSamples int, silent bool) {
	var records map[mutations.Mutation]Count
	var silence string
	if silent {
		records = c.possibleSilent
		silence = "silent"
	} else {
		records = c.possibleNS
		silence = "non-silent"
	}

	found, total := 0, 0
	for k, v := range records {
		// fmt.Printf("%c%d%c: %d\n", k.From, k.Pos+1, k.To, v)
		_, _, err := c.ref.Orfs.GetCodonOffset(k.Pos)
		if err != nil {
			continue
		}

		if v.count >= minSamples {
			found++
		}
		total++
	}
	fmt.Printf("%d/%d possible %s (%.4f) were actually "+
		"found above depth %d\n",
		found, total, silence,
		float64(found)/float64(total), c.minDepth)
}

func (c *CountAll) DisplayEvolution() {
	type Result struct {
		date int
		mean float64
	}
	results := make([]Result, 0, len(c.evolution))
	for k, v := range c.evolution {
		results = append(results, Result{k, v.Get()})
	}
	utils.SortByKey(results, false, func(r Result) int {
		return r.date
	})

	fd, fp := utils.WriteFile("evolution.dat")
	defer fd.Close()

	for _, r := range results {
		fmt.Fprintln(fp, r.date, r.mean)
	}
	fp.Flush()

	fmt.Println("Wrote evolution.dat")
}

func (c *CountAll) GraphPossible(silent bool, minSamples int, fname string) {
	var records map[mutations.Mutation]Count
	if silent {
		records = c.possibleSilent
	} else {
		records = c.possibleNS
	}

	fd, fp := utils.WriteFile(fname)
	defer fd.Close()

	type datum struct {
		pos           int
		silent        bool
		count         int
		maxDepthRatio float64
	}
	data := make([]datum, 0, len(records))
	for k, v := range records {
		if v.count >= minSamples {
			data = append(data,
				datum{k.Pos, k.Silent, v.count, v.maxDepthRatio})
		}
	}
	utils.SortByKey(data, true, func(d datum) int {
		return d.pos
	})

	for _, datum := range data {
		// fmt.Fprintln(fp, datum.pos, datum.count)
		fmt.Fprintln(fp, datum.pos, datum.maxDepthRatio)
	}
	fp.Flush()
	fmt.Printf("Wrote %s\n", fname)
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

var ROOT string = "/fs/bowser/genomes/raw_reads/"

func ProcessReads(db *database.Database,
	ids []database.Id, prefix string,
	handler ReadHandler) {

	for _, id := range ids {
		record := db.Get(id)
		pu := FindPileup(record, path.Join(ROOT, prefix))
		if pu == nil {
			continue
		}
		handler.Process(record, pu)
	}
}

func ProcessAll(dir string, handler ReadHandler) {
	matches, _ := filepath.Glob(path.Join(ROOT, dir, "*.txt.gz"))
	for _, m := range matches {
		pu, err := pileup.Parse2(m)
		if err != nil {
			log.Fatal(err)
		}
		handler.Process(nil, pu)
	}
}

func DisplayReads(db *database.Database,
	ids []database.Id, minDepth int, prefix string) {
	fmt.Println("Date AccNo SRA Region Country class 8782:C|T 28144:C|T QSlocs")

	var d Display
	d.Init(minDepth)
	ProcessReads(db, ids, prefix, &d)
}

func (c *CountAll) AverageMDR() float64 {
	var total, count float64
	for _, v := range c.counts {
		total += v.maxDepthRatio
		count++
	}
	return total / count
}

func CountEverything(db *database.Database,
	class string, dir string,
	ids []database.Id, minDepth int, minSamples int,
	prefix string, cutoff time.Time, nonMaj bool, transitions bool) {
	var c CountAll
	ref := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	c.Init(ref, minDepth, cutoff, nonMaj)

	if class == "all" {
		ProcessAll(dir, &c)
	} else {
		ProcessReads(db, ids, dir, &c)
	}

	c.SortDates()
	c.Display(minSamples)
	c.DisplayPossible(minSamples, false)
	c.DisplayPossible(minSamples, true)
	fmt.Printf("Average max depth ratio: %.2f\n", c.AverageMDR())

	c.DisplayEvolution()

	c.GraphPossible(false, minSamples, class+"-NS.dat")
	c.GraphPossible(true, minSamples, class+"-S.dat")

	if transitions {
		c.TransitionRatios()
	}

	gpiName := class + "-plot.gpi"
	fd, fp := utils.WriteFile(gpiName)
	defer fd.Close()
	fmt.Fprintf(fp,
		"set title \"%s\"\n", class)
	fmt.Fprintf(fp,
		"plot \"%s-NS.dat\" with impulses title \"NS\", "+
			"\"%s-S.dat\" with impulses title \"S\"\n",
		class, class)
	fp.Flush()
	fmt.Printf("Wrote %s\n", gpiName)
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
		minSamples    int
		findSequences bool
		class         string
		count         bool
		countCT       bool
		cutoffS       string
		cutoff        time.Time
		locations     bool
		showPoss      bool
		nonMaj        bool
		transitions   bool
		dir           string
	)

	flag.BoolVar(&findSequences, "f", false, "Find the sequences")
	flag.StringVar(&class, "class", "first100", "Classes to show")
	flag.IntVar(&minDepth, "min-depth", 3, "Min depth")
	flag.IntVar(&minSamples, "min-samp", 5, "Min samples")
	flag.BoolVar(&count, "countCT", false, "Count silent CT")
	flag.BoolVar(&count, "count", false, "Count everything")
	flag.StringVar(&cutoffS, "cutoff", "", "Count before date")
	flag.BoolVar(&locations, "locations", false, "Just show locations")
	flag.BoolVar(&showPoss, "show-poss", false, "Show possible silent")
	flag.BoolVar(&nonMaj, "non-maj", false, "Only show non-majority")
	flag.BoolVar(&transitions, "show-trans", false, "Show transitions")
	flag.StringVar(&dir, "dir", "pileups", "Directory")

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
		CountEverything(db, class, dir,
			records, minDepth, minSamples, class, cutoff, nonMaj, transitions)
	} else {
		DisplayReads(db, records, minDepth, class)
	}
}
