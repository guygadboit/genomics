package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/pileup"
	"genomics/stats"
	"log"
	"slices"
	"strings"
)

func printSorted(counts map[int]int) {
	type record struct {
		k, v int
	}
	records := make([]record, 0, len(counts))
	total := 0
	for k, v := range counts {
		records = append(records, record{k, v})
		total += v
	}
	slices.SortFunc(records, func(a, b record) int {
		if a.v < b.v {
			return 1
		}
		if a.v > b.v {
			return -1
		}
		return 0
	})

	fmt.Printf("Total matches: %d\n", total)
	for _, rec := range records {
		fmt.Printf("%d: %d\n", rec.k, rec.v)
	}
}

func majority(alleles map[byte]int) byte {
	var ret byte
	best := -1
	for k, v := range alleles {
		if v > best {
			best = v
			ret = k
		}
	}
	return ret
}

func iterateSilent(g *genomes.Genomes, cb func(int, byte)) {
	for _, mut := range mutations.PossibleSilentMuts(g, 0) {
		cb(mut.Pos, mut.To)
	}
}

func iterateAll(g *genomes.Genomes, cb func(int, byte)) {
	for i := 0; i < g.Length(); i++ {
		for _, nt := range []byte{'G', 'A', 'T', 'C'} {
			if nt != g.Nts[0][i] {
				cb(i, nt)
			}
		}
	}
}


func ExpectedMajorityRate(g *genomes.Genomes,
	requireSilent bool) (int, int, float64) {
	var iterate func(g *genomes.Genomes, cb func(int, byte))
	if requireSilent {
		iterate = iterateSilent
	} else {
		iterate = iterateAll
	}

	var count, total int
	iterate(g, func(pos int, nt byte) {
		alleles := make(map[byte]int)
		for i := 1; i < g.NumGenomes(); i++ {
			alleles[g.Nts[i][pos]]++
		}
		if majority(alleles) == nt {
			count++
		}
		total++
	})

	fmt.Println(count, total)
	return count, total-count, float64(count)/float64(total)
}

func Compare(pileup *pileup.Pileup,
	g *genomes.Genomes, minDepth int, requireSilent bool) {
	counts := make(map[int]int)

	// The total number of differences from g.Nts[0] with at least minDepth,
	// and also silent if requireSilent.
	diffs := 0

	// The number of diffs that match the majority allele in the outgroup
	totalMaj := 0

	for i := 0; i < g.Length(); i++ {
		rec := pileup.Get(i)
		if rec == nil {
			continue
		}
		for rank, read := range rec.Reads {
			if read.Nt == g.Nts[0][i] {
				continue
			}
			if read.Depth <= minDepth {
				continue
			}
			silent, _, _ := genomes.IsSilentWithReplacement(g,
				i, 0, 0, []byte{read.Nt})
			if requireSilent && !silent {
				continue
			}

			diffs++
			matches := make([]string, 0)
			alleles := make(map[byte]int)

			for j := 1; j < g.NumGenomes(); j++ {
				alleles[g.Nts[j][i]]++
				if read.Nt == g.Nts[j][i] {
					matches = append(matches, fmt.Sprintf("%d", j))
					counts[j]++
				}
			}

			var silentS string
			if silent {
				silentS = "*"
			}

			var majS string
			if read.Nt == majority(alleles) {
				majS = "M"
				totalMaj++
			}

			fmt.Printf("%c%d%c%s%s depth:%d rank:%d matches:%d:%s ",
				g.Nts[0][i], rec.Pos+1, read.Nt, silentS, majS, read.Depth,
				rank, len(matches), strings.Join(matches, ","))

			for k, v := range alleles {
				fmt.Printf("%c:%d ", k, v)
			}
			fmt.Printf("\n")
		}
	}
	rate := float64(totalMaj)/float64(diffs)
	printSorted(counts)

	a, b, rate := totalMaj, diffs - totalMaj, float64(totalMaj)/float64(diffs)
	c, d, expectedRate := ExpectedMajorityRate(g, requireSilent)

	fmt.Printf("%d/%d %.2f are majority\n", totalMaj, diffs, rate)
	fmt.Printf("Expected majority rate: %.2f\n", expectedRate)

	var ct stats.ContingencyTable
	ct.Init(a, b, c, d)
	OR, p := ct.FisherExact(stats.GREATER)
	fmt.Printf("OR=%.2f p=%g\n", OR, p)
}

func main() {
	var (
		fasta, orfs string
		minDepth    int
		silent      bool
	)

	flag.StringVar(&fasta, "fasta", "", "Reference alignment")
	flag.StringVar(&orfs, "orfs", "", "Reference ORFs")
	flag.IntVar(&minDepth, "min-depth", 4, "Minimum depth")
	flag.BoolVar(&silent, "silent", false, "Require silent")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	for _, arg := range flag.Args() {
		pileup, err := pileup.Parse(arg)
		if err != nil {
			log.Fatal(err)
		}
		Compare(pileup, g, minDepth, silent)
	}
}
