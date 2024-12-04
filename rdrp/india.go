package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
)

// These are those Indian sequences with the 42 NS muts in the RdRP
var interesting []string = []string{
	"EPI_ISL_10716753",
	"EPI_ISL_10716749",
	"EPI_ISL_10716805",
	"EPI_ISL_10716754",
	"EPI_ISL_10716752",
	"EPI_ISL_10716750",
}

// Get the intersection of the AA Changes in these sequences
func GetAAChanges(db *database.Database) map[database.AAMutation]bool {
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
	return muts
}

func ShowSequences(db *database.Database) {
	ids := db.GetByAccession(interesting...)

	muts := GetAAChanges(db)

	matches := utils.FromSet(db.SearchByAAMut(utils.FromSet(muts), len(muts)))
	db.Sort(matches, database.COLLECTION_DATE)

	for _, match := range matches {
		r := db.Get(match)
		fmt.Println(r.ToString(), r.Host)
	}

	// OK now look for how many other sequences in the DB share anything with
	// these guys, sorted by the most matches first. YOU ARE HERE

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
