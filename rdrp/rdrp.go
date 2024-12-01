package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
	"strings"
	"slices"
)

type RdRPMutation struct {
	database.Mutation
	ids	[]database.Id
}

func RdRPVariants(db *database.Database) {
	// The RdRP is nsp7, nsp8 and nsp12. Find all sequences with variations
	// anywhere in there.
	positions := make([]utils.OneBasedPos, 0)

	addRange := func(start, end utils.OneBasedPos) {
		for i := start; i <= end; i++ {
			positions = append(positions, i)
		}
	}

	addRange(11843, 12091)	// nsp7
	addRange(12094, 12685)	// nsp8
	addRange(13442, 16236)	// nsp12

	posSet := utils.ToSet(positions)
	ids := db.SearchByMutPosition(positions, 1)

	// Now group them by the actual unique mutations. Make an index of mutation
	// to ids that have it.
	index := make(map[database.Mutation]database.IdSet)

	for id, _ := range ids {
		r := db.Get(id)
		for _, mut := range r.NucleotideChanges {

			// Don't print out anything referring to other muts these records
			// may have that aren't in nsp7/8/12
			if !posSet[mut.Pos] {
				continue
			}
			_, there := index[mut]
			if !there {
				index[mut] = make(database.IdSet)
			}
			index[mut][id] = true
		}
	}

	// Now make that into something we can sort
	muts := make([]RdRPMutation, 0, len(index))
	for k, v := range index {
		muts = append(muts, RdRPMutation{k, utils.FromSet(v)})
	}

	slices.SortFunc(muts, func(a, b RdRPMutation) int {
		ka, kb := len(a.ids), len(b.ids)
		if ka < kb {
			return 1
		}
		if ka > kb {
			return -1
		}
		return 0
	})

	for _, mut := range muts {
		fmt.Printf("%s in %d sequences\n", mut.ToString(), len(mut.ids))

		db.Sort(mut.ids, database.COLLECTION_DATE)
		for _, id := range mut.ids {
			r := db.Get(id)
			fmt.Println(r.ToString())
		}
	}
}

func main() {
	db := database.NewDatabase()
	RdRPVariants(db)
}
