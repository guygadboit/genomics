package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
	"time"
)

func main() {
	db := database.NewDatabase()

	muts := database.ParseMutations("T18060C,T8782C,C28144T")	// α1, α2, α3
	// muts := database.ParseMutations("A17858G,C17747T")	// ν1 and ν2
	ids := db.SearchByMuts(muts, len(muts))
	fmt.Printf("Found %d with those muts\n", len(ids))

	/*
	cutoff := utils.Date(2020, 1, 30)
	ids = db.Filter(ids, func(r *database.Record) bool {
		return r.CollectionDate.Compare(cutoff) < 0
	})
	*/

	sorted := utils.FromSet(ids)
	db.Sort(sorted, database.COLLECTION_DATE)

	for _, id := range sorted {
		r := &db.Records[id]
		// fmt.Println(r.Summary())
		fmt.Println(r.CollectionDate.Format(time.DateOnly),
			r.SubmissionDate.Format(time.DateOnly),
			r.GisaidAccession, r.Country, r.Region)
	}
}
