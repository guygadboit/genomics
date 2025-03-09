package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
	"log"
	"strings"
	"slices"
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

func Compare(pileup *pileup.Pileup, g *genomes.Genomes, minDepth int) {
	counts := make(map[int]int)

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

			matches := make([]string, 0)

			for j := 1; j < g.NumGenomes(); j++ {
				if read.Nt == g.Nts[j][i] {
					matches = append(matches, fmt.Sprintf("%d", j))
					counts[j]++
				}
			}

			if len(matches) != 0 {
				fmt.Printf("%d%c depth:%d rank:%d matches:%s\n", rec.Pos+1,
					read.Nt, read.Depth, rank, strings.Join(matches, ","))
			}
		}
	}
	printSorted(counts)
}

func main() {
	var (
		fasta, orfs string
		minDepth    int
	)

	flag.StringVar(&fasta, "fasta", "", "Reference alignment")
	flag.StringVar(&orfs, "orfs", "", "Reference ORFs")
	flag.IntVar(&minDepth, "min-depth", 4, "Minimum depth")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	for _, arg := range flag.Args() {
		pileup, err := pileup.Parse(arg)
		if err != nil {
			log.Fatal(err)
		}
		Compare(pileup, g, minDepth)
	}
}
