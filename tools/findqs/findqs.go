package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
	"log"
	"strings"
)

func Compare(pileup *pileup.Pileup, g *genomes.Genomes, minDepth int) {
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
				}
			}

			// Count up which others were matched most often YOU ARE HERE

			if len(matches) != 0 {
				fmt.Printf("%d%c depth:%d rank:%d matches:%s\n", rec.Pos+1,
					read.Nt, read.Depth, rank, strings.Join(matches, ","))
			}
		}
	}
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
