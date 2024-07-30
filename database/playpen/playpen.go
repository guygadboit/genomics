package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
)

func main() {
	db := database.NewDatabase()

	ids := db.Filter(nil, func(r *database.Record) bool {
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
		fmt.Println(r.Summary())
	}
}
