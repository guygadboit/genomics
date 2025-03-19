package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"slices"
)

func SortedSimilarity(g *genomes.Genomes) {
	type result struct {
		index int
		ss    float64
	}

	results := make([]result, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		results[i] = result{i, g.SequenceSimilarity(i, 0)}
	}

	slices.SortFunc(results, func(a, b result) int {
		if a.ss < b.ss {
			return 1
		}
		if a.ss > b.ss {
			return -1
		}
		return 0
	})

	for _, r := range results {
		fmt.Printf("%d: %s %.2f%%\n", r.index, g.Names[r.index], r.ss*100)
	}
}

func ExhaustiveSimilarity(g *genomes.Genomes) {
	type result struct {
		a, b int
		ss   float64
	}

	results := make([]result, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		for j := i + 1; j < g.NumGenomes(); j++ {
			results = append(results, result{i, j, g.SequenceSimilarity(i, j)})
		}
	}

	slices.SortFunc(results, func(a, b result) int {
		if a.ss < b.ss {
			return 1
		}
		if a.ss > b.ss {
			return -1
		}
		return 0
	})

	for _, r := range results {
		fmt.Printf("%d %s,%d %s: %.2f%%\n",
			r.a, g.Names[r.a], r.b, g.Names[r.b], r.ss*100)
	}
}

func main() {
	var include, exclude, outName string
	var summary, ss, sss, ess bool
	var removeGaps bool

	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.StringVar(&exclude, "e", "", "Genomes to exclude")
	flag.BoolVar(&summary, "s", false, "Summary")
	flag.BoolVar(&ss, "ss", false, "Similiarity summary")
	flag.BoolVar(&sss, "sss", false, "Sorted similiarity summary")
	flag.BoolVar(&ess, "ess", false, "Exhaustive similarity summary")
	flag.BoolVar(&removeGaps, "g", false, "Remove Gaps")
	flag.StringVar(&outName, "o", "filtered.fasta", "Output filename")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", false)

	if summary {
		fmt.Printf("%d nucleotides\n", g.Length())
		g.PrintSummary()
		return
	}

	if ss {
		for i, n := range g.Names {
			ss := g.SequenceSimilarity(i, 0)
			fmt.Printf("%d: %s %.2f%%\n", i, n, ss*100)
		}
		return
	} else if sss {
		SortedSimilarity(g)
		return
	} else if ess {
		ExhaustiveSimilarity(g)
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

	var g2 *genomes.Genomes
	if len(which) > 0 {
		g2 = g.Filter(which...)
	} else {
		g2 = g
	}
	if removeGaps {
		g2.RemoveGaps()
	}
	g2.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
