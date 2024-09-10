package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
	"time"
)

func main() {
	db := database.NewDatabase()

	cutoff := utils.Date(2020, 1, 30)
	ids := db.Filter(nil, func(r *database.Record) bool {
		return r.CollectionDate.Compare(cutoff) < 0

		if r.Host != "Human" {
			return false
		}
		if r.Country != "Egypt" {
			return false
		}

		for _, mut := range r.AAChanges {
			if mut.Gene == "ORF1a" && (mut.Pos == 124 || mut.Pos == 125) {
				fmt.Println(mut)
				return true
			}
		}

		return false
	})

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
