package main

import (
	"fmt"
	"genomics/stats"
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
		results := stats.Blast(bc, genome, ins.Nts, 1, 1, stats.NOT_VERBOSE)
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
