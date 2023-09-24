package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"io"
	"log"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"
)

type Insertion struct {
	id      int    // The order we found them in. Should be a line number
	pos     int    // Where
	nts     []byte // What
	nSeqs   int    // How many times
	inWH1   bool   // Is it in WH1?
	inHuman bool   // Is it in human? Note that most short things will be
}

func (i *Insertion) ToString() string {
	return fmt.Sprintf("%s at %d (%d seqs)", string(i.nts),
		i.pos, i.nSeqs)
}

func atoi(s string) int {
	ret, err := strconv.Atoi(s)
	if err != nil {
		log.Fatal("Bad integer")
	}
	return ret
}

func LoadInsertions(fname string, minLen int, minSeqs int) []Insertion {
	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal("Can't open %s", fname)
	}
	defer fd.Close()

	ret := make([]Insertion, 0)
	fp := bufio.NewReader(fd)

	// Match for [GATC]+, so ignore any with Ns or weird ambiguous nts like HDK
	// etc.
	pat := regexp.MustCompile(`ins_(\d+):([GATC]+) \((\d+) seqs\)`)
	id := 1

reading:
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break reading
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		id++
		line = strings.TrimSpace(line)
		groups := pat.FindAllStringSubmatch(line, -1)
		if groups == nil {
			continue
		}

		ins := Insertion{id,
			atoi(groups[0][1]),
			[]byte(groups[0][2]),
			atoi(groups[0][3]),
			false, false}

		if len(ins.nts) < minLen {
			continue
		}

		if ins.nSeqs < minSeqs {
			continue
		}

		// Some of the insertions seem to be nearly all A. This looks bogus. So
		// filter them out.
		var aCount int
		for i := 0; i < len(ins.nts); i++ {
			if ins.nts[i] == 'A' {
				aCount++
			}
		}
		if (float64(aCount) / float64(len(ins.nts))) > 0.99 {
			continue
		}

		ret = append(ret, ins)
	}

	return ret
}

func Summary(insertions []Insertion) {
	fmt.Printf("%d insertions\n", len(insertions))

	var total int
	for i := 0; i < len(insertions); i++ {
		total += len(insertions[i].nts)
	}
	fmt.Printf("%d nts altogether (average length %.2f)\n",
		total, float64(total)/float64(len(insertions)))
}

/*
	Call cb for every insertion found in nts forwards or backwards
*/
func search(insertion *Insertion, genome *genomes.Genomes, which int,
	tol float64, cb func(*Insertion)) {

	var search genomes.Search
	nts := insertion.nts

	for search.Init(genome, 0, nts, tol); !search.End(); search.Next() {
		cb(insertion)
		return
	}

	rc := utils.ReverseComplement(insertion.nts)
	for search.Init(genome, 0, rc, tol); !search.End(); search.Next() {
		cb(insertion)
		return
	}
}

func findInVirus(name string,
	insertions []Insertion, minLength int, mark bool, tol float64) {
	virus := genomes.LoadGenomes(fmt.Sprintf("../fasta/%s.fasta", name),
		"", false)

	fname := fmt.Sprintf("%s-insertions.txt", name)
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	reportFound := func(ins *Insertion) {
		// fmt.Fprintf(w, "ins_%d:%s (%d seqs)\n", ins.pos, ins.nts, ins.nSeqs)
		if mark {
			ins.inWH1 = true
		}
	}

	var found, count int

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]

		nts := ins.nts
		if len(nts) < minLength {
			continue
		}
		count++

		search(ins, virus, 0, tol, func(ins *Insertion) {
			reportFound(ins)
			found++
		})
	}

	w.Flush()
	fmt.Printf("Length %d: %d (/%d) were found in %s\n", minLength,
		found, count, name)
}

/*
	Mark the ones you find in human, only if they aren't already in WH1
*/
func findInHuman(insertions []Insertion, minLength int, tol float64) {
	fmt.Printf("Loading...\n")
	g := genomes.LoadGenomes("/fs/f/genomes/human/GRCh38_latest_genomic.fna.gz",
		"", true)
	fmt.Printf("Loaded human\n")

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		nts := ins.nts

		if len(nts) < minLength {
			continue
		}

		if ins.inWH1 {
			continue
		}

		search(ins, g, 0, tol, func(ins *Insertion) {
			fmt.Printf("%d (length %d) is in human\n", ins.id, len(ins.nts))
			ins.inHuman = true
		})

		if i%10 == 0 {
			fmt.Printf("Searched %d/%d\n", i, len(insertions))
		}
	}
}

/*
	Ponderous sort but less trouble than updating my whole system to the latest
	Go version that has slices.SortFunc
*/

type posCount struct {
	pos   int
	count int
}

type posCounts []posCount

func (p posCounts) Len() int {
	return len(p)
}

func (p posCounts) Less(i, j int) bool {
	return p[i].pos < p[j].pos
}

func (p posCounts) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func byLocation(insertions []Insertion, minLength int) {
	positions := make(map[int]int)

	for i := 0; i < len(insertions); i++ {
		count, _ := positions[insertions[i].pos]
		positions[insertions[i].pos] = count + 1
	}

	posCounts := make(posCounts, 0, len(positions))
	for k, v := range positions {
		posCounts = append(posCounts, posCount{k, v})
	}

	sort.Sort(posCounts)

	fd, err := os.Create("locations.txt")
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	for i := 0; i < posCounts.Len(); i++ {
		pc := posCounts[i]
		fmt.Fprintf(w, "%d %d\n", pc.pos, pc.count)
	}

	fmt.Printf("Wrote locations.txt\n")
}

type filter int

const (
	ANYTHING filter = iota
	WH1_ONLY
	NOT_WH1_ONLY
)

/*
	Call cb on all the insertions that match the filter. This requires that
	you've already marked which ones are in WH1 with findInVirus.
*/
func filterInsertions(insertions []Insertion,
	minLength int, filter filter, cb func(*Insertion), verbose bool) {
	var name string

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]

		if len(ins.nts) < minLength {
			continue
		}

		switch filter {
		case WH1_ONLY:
			name = "In WH1"
			if !ins.inWH1 {
				continue
			}
		case NOT_WH1_ONLY:
			name = "Not in WH1"
			if ins.inWH1 {
				continue
			}
		case ANYTHING:
			name = "Anything"
			break
		}

		if verbose {
			fmt.Printf("%s: %d ins_%d %s (%d seqs)\n", name,
				len(ins.nts), ins.pos, string(ins.nts), ins.nSeqs)
		}

		cb(ins)
	}
}

/*
	Output them all in one fasta file with NNN between each of them. Useful if
	you want to look at dinucleotide composition etc.
*/
func outputCombinedFasta(fname string, name string,
	insertions []Insertion, minLength int, filter filter, verbose bool) int {
	sep := []byte("NNN")
	nts := make([]byte, 0)
	var count int

	cb := func(ins *Insertion) {

		nts = append(nts, ins.nts...)
		nts = append(nts, sep...)
		count++
	}

	filterInsertions(insertions, minLength, filter, cb, verbose)

	var orfs genomes.Orfs
	genomes := genomes.NewGenomes(orfs, 1)
	genomes.Nts[0] = nts
	genomes.Names[0] = name
	genomes.Save("SC2 Insertions", fname, 0)

	fmt.Printf("%s: %d insertions\n", name, count)
	return count
}

/*
	Put all the insertions matching the filters into a single fasta file with
	headers between each one. I think you should then be able to BLAST the
	whole lot.
*/
func outputFasta(fname string, name string,
	insertions []Insertion, minLength int, filter filter, verbose bool) {

	var orfs genomes.Orfs
	genomes := genomes.NewGenomes(orfs, 0)

	cb := func(ins *Insertion) {
		genomes.Nts = append(genomes.Nts, ins.nts)
		name := fmt.Sprintf("ins_%d_%d", ins.pos, ins.id)
		genomes.Names = append(genomes.Names, name)
	}

	filterInsertions(insertions, minLength, filter, cb, verbose)
	genomes.SaveMulti(fname)
}

func outputDinucs(fname string,
	insertions []Insertion, minLength int, filter filter, verbose bool) {

	fd, err := os.Create(fname)
	if err != nil {
		log.Fatalf("Can't create %s", fname)
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)
	fmt.Fprintf(w, "id len G+C, CpG TpA insert\n")

	cb := func(ins *Insertion) {
		dp := genomes.CalcProfile(ins.nts)
		fmt.Fprintf(w, "%d %d %.3f %.3f %.3f %s %t\n",
			ins.id, len(ins.nts), dp.GC, dp.CpG, dp.TpA, string(ins.nts),
			ins.inHuman)
	}
	filterInsertions(insertions, minLength, filter, cb, verbose)
	w.Flush()
}

func showLength(insertions []Insertion) {
	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		if ins.inWH1 {
			continue
		}
		fmt.Printf("%d ins_%d %s\n", len(ins.nts), ins.pos, ins.nts)
	}
}

func otherHCoVs(insertions []Insertion, tol float64) {
	covs := [...]string{"229E", "NL63", "OC43", "HKU1", "WH1"}

	var wg sync.WaitGroup
	for i := 0; i < len(covs); i++ {
		wg.Add(1)
		go func(i int) {
			findInVirus(covs[i], insertions, 10, false, tol)
			wg.Done()
		}(i)
	}
	wg.Wait()
}

func main() {
	var tol float64
	var verbose bool

	flag.Float64Var(&tol, "tol", 0.0, "Tolerance")
	flag.BoolVar(&verbose, "v", false, "Verbose")
	flag.Parse()

	insertions := LoadInsertions("insertions.txt", 6, 2)
	findInVirus("WH1", insertions, 6, true, tol)
	findInHuman(insertions, 20, tol)

	utils.Sort(len(insertions), true,
		func(i, j int) {
			insertions[i], insertions[j] = insertions[j], insertions[i]
		},
		nil,
		func(i int) float64 {
			return float64(len(insertions[i].nts))
		})

	// showLength(insertions)
	// byLocation(insertions, 9)

	outputFasta("InsertionsNotFromWH1.fasta", "FromWH1",
		insertions, 6, NOT_WH1_ONLY, verbose)

	outputDinucs("InsertionsNotFromWH1.txt",
		insertions, 6, NOT_WH1_ONLY, verbose)
}
