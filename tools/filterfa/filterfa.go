package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
)

func main() {
	var include, exclude, outName string
	var summary, ss bool
	var removeGaps bool

	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.StringVar(&exclude, "e", "", "Genomes to exclude")
	flag.BoolVar(&summary, "s", false, "Summary")
	flag.BoolVar(&ss, "ss", false, "Similiarity summary")
	flag.BoolVar(&removeGaps, "g", false, "Remove Gaps")
	flag.StringVar(&outName, "o", "filtered.fasta", "Output filename")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", false)

	if summary {
		g.PrintSummary()
		return
	}

	if ss {
		for i, n := range g.Names {
			ss := g.SequenceSimilarity(i, 0)
			fmt.Printf("%d: %s %.2f%%\n", i, n, ss*100)
		}
		return
	}

	if include != "" && exclude != "" {
		log.Fatalf("Specify either include or exclude, not both")
	}

	var which []int

	if include != "" {
		which = utils.ParseInts(include, ",")
	}

	if exclude != "" {
		which = make([]int, 0)
		exc := utils.ToSet(utils.ParseInts(exclude, ","))
		for i := 0; i < g.NumGenomes(); i++ {
			if !exc[i] {
				which = append(which, i)
			}
		}
	}

	g2 := g.Filter(which...)
	if removeGaps {
		g2.RemoveGaps()
	}
	g2.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
