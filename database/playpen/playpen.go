package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
)

func main() {
	db := database.NewDatabase()

	ids := db.SearchByMuts(database.ParseMutations("G11083T,G26144T"), 1)
	ids = db.Filter(ids, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		cutoff := utils.Date(2020, 1, 29)
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
