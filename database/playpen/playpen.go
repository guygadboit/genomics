package main

import (
	"fmt"
	"time"
	"genomics/database"
	"genomics/utils"
)

func main() {
	db := database.NewDatabase()

	// ids := db.SearchByMuts(database.ParseMutations("T5959G,A8651C,G16206A"), 1)
	ids := db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		if r.Country != "Egypt" {
			return false
		}
		return true
	})

	sorted := utils.FromSet(ids)
	db.Sort(sorted, database.COLLECTION_DATE)

	total, count := 0, 0
	for _, id := range sorted {
		r := &db.Records[id]
		fmt.Println(r.GisaidAccession,
			r.CollectionDate.Format(time.DateOnly), len(r.NucleotideChanges))
		total += len(r.NucleotideChanges)
		count++
	}
	fmt.Printf("Average %f\n", float64(total)/float64(count))
}
