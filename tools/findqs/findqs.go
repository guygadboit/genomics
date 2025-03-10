package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
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
	fmt.Printf("%d/%d %.2f are majority\n", totalMaj, diffs, rate)
	printSorted(counts)
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
