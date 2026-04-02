package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/utils"
)

func main() {
	var (
		ntMuts string
	)
	flag.StringVar(&ntMuts, "nt", "", "Nucleotide muts")
	flag.Parse()

	db := database.NewDatabase()
	muts := database.ParseMutations(ntMuts)

	positions := make([]utils.OneBasedPos, len(muts))
	for i, mut := range muts {
		positions[i] = mut.Pos
	}

	matches := make([]*database.Record, 0)
	search := database.NewNtMutationIndexSearch(db, positions)
	counts := make(map[database.Id]int)

	for i, target := range muts {
		ids, _ := search.Get(i)
		for id, _ := range ids {
			record := db.Get(id)
			if record.Host != "Human" {
				continue
			}
			got := record.HasMuts(database.Mutations{target})
			if len(got) != 1 {
				continue
			}
			counts[id]++
			matches = append(matches, record)
		}
	}

	type result struct {
		id database.Id
		count int
	}
	results := make([]result, 0, len(counts))
	for k, v := range counts {
		results = append(results, result{k, v})
	}

	utils.SortByKey(results, true, func(r result) int {
		return r.count
	})

	for _, r := range results {
		rec := db.Get(r.id)
		fmt.Printf("%d/%d: %s\n", r.count, len(muts), rec.Summary())
	}
}
