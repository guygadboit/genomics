package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/utils"
)

func main() {
	db := database.NewDatabase()

	flag.Parse()

	// Given a list of EPI_ISL numbers print out the SRAs for any we can find
	for _, fname := range flag.Args() {
		utils.Lines(fname, func(line string, err error) bool {
			ids := db.GetByAccession(line)
			for _, id := range ids {
				record := db.Get(id)
				for _, sra := range record.SRA {
					fmt.Println(sra)
				}
			}
			return true
		})
	}
}
