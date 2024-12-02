package main

import (
	"fmt"
	"genomics/database"
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

func ShowSequences(db *database.Database) {
	ids := db.GetByAccession(interesting...)

	/*
		OK I guess show what they have in common, especially AA changes in the
		RdRP. Then we want to look for subsets of those elsewhere is the idea,
		including if they ever crop up but silently mutated. YOU ARE HERE.
	*/
	for _, id := range ids {
		r := db.Get(id)
		fmt.Println(r.GisaidAccession,
			len(r.AAChanges), len(r.NucleotideChanges))
	}
}
