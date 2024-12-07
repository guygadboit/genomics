package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
)

var interesting []string = []string{
	"EPI_ISL_1133157",
}

// Get the intersection of the AA Changes in these sequences. Return them
// indexed by gene.
func GetAAChanges(db *database.Database) map[string]database.AAMutations {
	ids := db.GetByAccession(interesting...)

	var muts map[database.AAMutation]bool
	for _, id := range ids {
		r := db.Get(id)
		if muts == nil {
			muts = utils.ToSet(r.AAChanges)
		} else {
			muts = utils.Intersection(muts, utils.ToSet(r.AAChanges))
		}
	}

	ret := make(map[string]database.AAMutations)
	for mut, _ := range muts {
		if _, there := ret[mut.Gene]; !there {
			ret[mut.Gene] = make(database.AAMutations, 0)
		}
		ret[mut.Gene] = append(ret[mut.Gene], mut)
	}
	return ret
}

func ShowSequences(db *database.Database) {
	ids := db.GetByAccession(interesting...)

	allMuts := GetAAChanges(db)
	for gene, muts := range allMuts {
		mutsSet := utils.ToSet(muts)
		fmt.Printf("Gene %s has %d muts\n", gene, len(muts))

		for _, mut := range muts {
			fmt.Printf("%s ", mut.ToString())
		}
		fmt.Println()

		matches := db.SearchByAAMut(muts, 1)
		for _, match := range matches {
			r := db.Get(match.Id)
			ourMuts := utils.ToSet(r.AAChanges)
			ix := utils.FromSet(utils.Intersection(mutsSet, ourMuts))

			if len(ix) == 0 {
				continue
			}

			fmt.Printf("%s (%s) %d\n", r.ToString(), r.Host, len(ix))

			utils.SortByKey(ix, func(mut database.AAMutation) int {
				return int(byte(mut.Gene[0]))*10000 + int(mut.Pos)
			}, false)
			for _, mut := range ix {
				fmt.Printf("%s ", mut.ToString())
			}
			fmt.Println()

		}
		fmt.Println()
	}

	return

	/*
		OK I guess show what they have in common, especially AA changes in the
		RdRP. Then we want to look for subsets of those elsewhere is the idea,
		including if they ever crop up but silently mutated.
	*/
	for _, id := range ids {
		r := db.Get(id)
		fmt.Println(r.GisaidAccession,
			len(r.AAChanges), len(r.NucleotideChanges))
	}
}
