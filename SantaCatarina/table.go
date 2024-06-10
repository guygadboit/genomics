package main

import (
	"bufio"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/utils"
	"slices"
	"time"
)

// Find all the muts present in any of the records here, as a sorted list of
// Mutations
func findMuts(db *database.Database, ids []database.Id) database.Mutations {
	muts := make(map[database.Mutation]bool)

	for _, id := range ids {
		for _, mut := range db.Records[id].NucleotideChanges {
			muts[mut] = true
		}
	}

	ret := make(database.Mutations, 0, len(muts))
	for k, _ := range muts {
		ret = append(ret, k)
	}
	slices.SortFunc(ret, func(a, b database.Mutation) int {
		return a.Pos - b.Pos
	})
	return ret
}

func allele(r *database.Record, mut database.Mutation) byte {
	if r == nil {
		return mut.From
	}
	has := r.HasMuts(database.Mutations{mut})
	if len(has) == 0 {
		return mut.From
	} else {
		return mut.To
	}
}

func insertFile(w *bufio.Writer, fname string) {
	utils.Lines(fname, func(line string, err error) bool {
		fmt.Fprintln(w, line)
		return true
	})
}

func insertMuts(muts database.Mutations,
	ins ...database.Mutation) database.Mutations {
	for _, mut := range ins {
		for i, m := range muts {
			if mut.Pos < m.Pos {
				muts = utils.Insert(muts, i, 1)
				muts[i] = mut
				break
			}
		}
	}
	return muts
}

// Make a table of who has what like the one in the Civet paper.
func MakeTable(db *database.Database,
	g *genomes.Genomes, ids []database.Id,
	muts database.Mutations,
	w *bufio.Writer) {
	insertFile(w, "head.html")

	if muts == nil {
		muts = findMuts(db, ids)
	}

	muts = insertMuts(muts,
		database.Mutation{8782, 'C', 'T', database.UNKNOWN},
		database.Mutation{28144, 'T', 'C', database.UNKNOWN},
		database.Mutation{18060, 'C', 'T', database.UNKNOWN},
	)

	fmt.Fprintf(w, "<tr id=\"headings\">")
	for i := 0; i < 5; i++ {
		fmt.Fprintf(w, "<th></th>\n")
	}
	for _, mut := range muts {
		fmt.Fprintf(w, "<th>%d</th>\n", mut.Pos)
	}
	fmt.Fprintf(w, "</tr>")

	insertFile(w, "florianopolis.html")

	reference := make([]byte, len(muts))
	outputMuts := func(ref bool, r *database.Record, muts database.Mutations) {
		for i, mut := range muts {
			nt := allele(r, mut)
			class := "same"
			if ref {
				reference[i] = nt
			} else {
				if reference[i] != nt {
					class = "different"
				}
			}
			fmt.Fprintf(w, "<td class=\"nt %s\">%c</td>\n", class, nt)
		}
	}

	for i, id := range ids {
		fmt.Fprintf(w, "<tr>")
		r := &db.Records[id]
		fmt.Fprintf(w, "<td>%s</td>", r.GisaidAccession)
		fmt.Fprintf(w, "<td>%s</td>", r.CollectionDate.Format(time.DateOnly))
		fmt.Fprintf(w, "<td>%s</td>", r.Region)
		fmt.Fprintf(w, "<td>%s</td>", r.City)
		fmt.Fprintf(w, "<td>%s</td>", r.Country)

		outputMuts(i == 0, r, muts)
		fmt.Fprintf(w, "</tr>")
	}

	fmt.Fprintf(w, "<tr><td>Wuhan-Hu-1</td>\n")
	for i := 0; i < 3; i++ {
		fmt.Fprintf(w, "<td></td>")
	}
	fmt.Fprintf(w, "<td>China</td>")

	outputMuts(false, nil, muts)
	fmt.Fprintf(w, "</tr>")
	insertFile(w, "tail.html")
}
