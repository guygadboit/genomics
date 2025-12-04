package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/utils"
	"time"
)

func main() {
	var withDates bool

	flag.BoolVar(&withDates, "d", false, "output the dates as well")
	flag.Parse()
	db := database.NewDatabase()

	// Given a list of EPI_ISL numbers print out the SRAs for any we can find
	for _, fname := range flag.Args() {
		utils.Lines(fname, func(line string, err error) bool {
			ids := db.GetByAccession(line)
			for _, id := range ids {
				record := db.Get(id)
				for _, sra := range record.SRA {
					if withDates {
						fmt.Printf("%s,%s\n", sra, record.CollectionDate.Format(
							time.DateOnly))
					} else {
						fmt.Println(sra)
					}
				}
			}
			return true
		})
	}
}
