package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
)

func main() {
	db := database.NewDatabase()

	ids := db.SearchByMuts(database.ParseMutations("T5959G,A8651C,G16206A"), 1)
	ids = db.Filter(ids, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		cutoff := utils.Date(2020, 3, 1)
		if r.CollectionDate.Compare(cutoff) > 0 {
			return false
		}
		return true
	})

	sorted := utils.FromSet(ids)
	db.Sort(sorted, database.COLLECTION_DATE)

	for _, id := range sorted {
		fmt.Println(db.Records[id].Summary())
	}
}
