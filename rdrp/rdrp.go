package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
	"slices"
)

type RdRPMutation struct {
	database.Mutation               // A mutation in the RdRP
	Ids               []database.Id // The ids that have it, sorted by date
}

type RdRPMutations struct {
	Db   *database.Database
	Muts []RdRPMutation

	// What RdRP muts does each id (that has any) have?
	ById map[database.Id][]database.Mutation
}

func NewRdRPMutations(db *database.Database) *RdRPMutations {
	// The RdRP is nsp7, nsp8 and nsp12. Find all sequences with variations
	// anywhere in there.
	positions := make([]utils.OneBasedPos, 0)

	addRange := func(start, end utils.OneBasedPos) {
		for i := start; i <= end; i++ {
			positions = append(positions, i)
		}
	}

	addRange(11843, 12091) // nsp7
	addRange(12094, 12685) // nsp8
	addRange(13442, 16236) // nsp12

	posSet := utils.ToSet(positions)
	matches := db.SearchByMutPosition(positions, 1)

	// Now group them by the actual unique mutations. Make an index of mutation
	// to ids that have it.
	index := make(map[database.Mutation]database.IdSet)

	for _, match := range matches {
		r := db.Get(match.Id)
		for _, mut := range r.NucleotideChanges {

			// Don't save anything referring to other muts these records may
			// have that aren't in nsp7/8/12
			if !posSet[mut.Pos] {
				continue
			}
			_, there := index[mut]
			if !there {
				index[mut] = make(database.IdSet)
			}
			index[mut][match.Id] = true
		}
	}

	// Now make that into something we can sort
	muts := make([]RdRPMutation, 0, len(index))
	for k, v := range index {
		muts = append(muts, RdRPMutation{k, utils.FromSet(v)})
	}

	// Sort them by the ones that occur the most often.
	slices.SortFunc(muts, func(a, b RdRPMutation) int {
		ka, kb := len(a.Ids), len(b.Ids)
		if ka < kb {
			return 1
		}
		if ka > kb {
			return -1
		}
		return 0
	})

	// And within each RdRPMutation, sort the ids by date
	for i, _ := range muts {
		db.Sort(muts[i].Ids, database.COLLECTION_DATE)
	}

	ret := &RdRPMutations{db, muts, nil}
	ret.findById()
	return ret
}

func (m *RdRPMutations) findById() {
	m.ById = make(map[database.Id][]database.Mutation)
	for _, rmut := range m.Muts {
		for _, id := range rmut.Ids {
			if _, there := m.ById[id]; !there {
				m.ById[id] = make([]database.Mutation, 0)
			}
			m.ById[id] = append(m.ById[id], rmut.Mutation)
		}
	}
}

func (m *RdRPMutations) Print() {
	for _, mut := range m.Muts {
		fmt.Printf("%s in %d sequences\n", mut.ToString(), len(mut.Ids))
		for _, id := range mut.Ids {
			r := m.Db.Get(id)
			fmt.Println(r.ToString())
		}
	}
}

func (m *RdRPMutations) PrintById() {
	type record struct {
		id   database.Id
		muts []database.Mutation
	}
	records := make([]record, 0)

	for k, v := range m.ById {
		nsMuts := make([]database.Mutation, 0)
		for _, mut := range v {
			if mut.Silence == utils.NON_SILENT {
				nsMuts = append(nsMuts, mut)
			}
		}
		records = append(records, record{k, nsMuts})
	}

	slices.SortFunc(records, func(a, b record) int {
		ka, kb := len(a.muts), len(b.muts)
		if ka < kb {
			return 1
		}
		if ka > kb {
			return -1
		}
		return 0
	})
	for _, rec := range records {
		r := m.Db.Get(rec.id)
		fmt.Printf("%s (%s): %d NS RdRP mutations\n",
			r.ToString(), r.Host, len(rec.muts))

		for _, mut := range rec.muts {
			fmt.Printf("%s ", mut.ToString())
		}
		fmt.Printf("\n")
	}
}

func main() {
	db := database.NewDatabase()
	ShowSequences(db)
	// muts := NewRdRPMutations(db)
	// muts.Print()
	// muts.PrintById()
}
