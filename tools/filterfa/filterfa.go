package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"slices"
)

func SortedSimilarity(g *genomes.Genomes, protein bool, ref int) {
	type result struct {
		index int
		ss    float64
	}

	results := make([]result, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		results[i] = result{i, g.SequenceSimilarity(i, ref, protein)}
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

func ExhaustiveSimilarity(g *genomes.Genomes, protein bool) {
	type result struct {
		a, b int
		ss   float64
	}

	results := make([]result, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		for j := i + 1; j < g.NumGenomes(); j++ {
			results = append(results,
				result{i, j, g.SequenceSimilarity(i, j, protein)})
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
	var (
		include, exclude, outName string
		summary, ss, sss, ess     bool
		removeGaps                bool
		protein                   bool
		dedupe                    bool
		densestFirst              bool
		ref                       int
	)

	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.StringVar(&exclude, "e", "", "Genomes to exclude")
	flag.BoolVar(&summary, "s", false, "Summary")
	flag.BoolVar(&ss, "ss", false, "Similiarity summary")
	flag.BoolVar(&sss, "sss", false, "Sorted similiarity summary")
	flag.BoolVar(&ess, "ess", false, "Exhaustive similarity summary")
	flag.BoolVar(&removeGaps, "g", false, "Remove Gaps")
	flag.BoolVar(&protein, "p", false, "Input is a protein")
	flag.StringVar(&outName, "o", "filtered.fasta", "Output filename")
	flag.BoolVar(&dedupe, "dd", false, "Remove duplicates")
	flag.BoolVar(&densestFirst, "df", false, "Put the densest first")
	flag.IntVar(&ref, "ref", 0, "Reference for similarity summary")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", false)

	if summary {
		fmt.Printf("%d nucleotides\n", g.Length())
		g.PrintSummary()
		return
	}

	if ss {
		for i, n := range g.Names {
			ss := g.SequenceSimilarity(i, ref, protein)
			fmt.Printf("%d: %s %.2f%%\n", i, n, ss*100)
		}
		return
	} else if sss {
		SortedSimilarity(g, protein, ref)
		return
	} else if ess {
		ExhaustiveSimilarity(g, protein)
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
	if dedupe {
		g2 = g.Dedupe()
	}
	if densestFirst {
		g2 = g.LeastGapsFirst()
	}

	g2.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
