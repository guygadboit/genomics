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
	var summary bool

	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.StringVar(&exclude, "e", "", "Genomes to exclude")
	flag.BoolVar(&summary, "s", false, "Summary")
	flag.StringVar(&outName, "o", "filtered.fasta", "Output filename")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", false)
	if summary {
		g.PrintSummary()
		return
	}

	if include != "" && exclude != "" {
		log.Fatalf("Specify either include or exclude, not both")
	}

	which := make([]int, 0)

	if include != "" {
		inc := utils.ToSet(utils.ParseInts(include))
		for i := 0; i < g.NumGenomes(); i++ {
			if inc[i] {
				which = append(which, i)
			}
		}
	}

	if exclude != "" {
		exc := utils.ToSet(utils.ParseInts(exclude))
		for i := 0; i < g.NumGenomes(); i++ {
			if !exc[i] {
				which = append(which, i)
			}
		}
	}

	g2 := g.Filter(which...)
	g2.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
