package main

import (
	"fmt"
	"genomics/stats"
	"genomics/utils"
	"log"
	"slices"
)

type BlastResult struct {
	stats.BlastResult
	id int
}

// Maps insertion id to a BlastResult. We will do one of these for each
// organism we're interested in.
func BlastInsertions(insertions []Insertion, genome string) []BlastResult {
	bc := stats.BlastDefaultConfig()
	ret := make([]BlastResult, 0)

	for _, ins := range insertions {
		results, _ := stats.Blast(bc, genome,
			ins.Nts, 1, 1, "", stats.NOT_VERBOSE)
		switch len(results) {
		case 0:
			continue
		case 1:
			fmt.Println(ins.Id, results[0].E)
			ret = append(ret, BlastResult{results[0], ins.Id})
		default:
			log.Fatal("Didn't think this should happen")
		}
	}

	slices.SortFunc(ret, func(a, b BlastResult) int {
		if a.E < b.E {
			return -1
		}
		if a.E > b.E {
			return 1
		}
		return 0
	})

	return ret
}

func ProkBlast(insertions []Insertion, filters []filterFunc) []BlastResult {
	ret := make([]BlastResult, 0)
	bc := stats.BlastDefaultConfig()
	filterInsertions(insertions, filters, func(ins *Insertion) {
		fmt.Printf("%s: %dnts\n", ins.ToString(), len(ins.Nts))
		results, err := stats.Blast(bc,
			"/fs/f/genomes/blast/prok/ref_prok_rep_genomes",
			ins.Nts, 10, 10, "", stats.NOT_VERBOSE)
		if err != nil {
			log.Print(err)
			return
		}
		for _, r := range results {
			fmt.Println(ins.Id, r.E, r.Organism)
			ret = append(ret, BlastResult{r, ins.Id})
		}

	}, false)
	return ret
}

func showResults(results []BlastResult) {
	organisms := make(map[string]int)
	for _, r := range results {
		organisms[r.Organism]++
		fmt.Printf("%d %s %d %f\n", r.id, r.Organism, r.Length, r.E)
	}

	type count struct {
		key   string
		value int
	}
	counts := make([]count, 0, len(results))

	for k, v := range organisms {
		counts = append(counts, count{k, v})
	}
	utils.SortByKey(counts, true, func(c count) int {
		return c.value
	})
	for _, c := range counts {
		fmt.Printf("%s: %d\n", c.key, c.value)
	}
}
