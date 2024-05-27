package main

import (
	"fmt"
	"genomics/database"
)

func main() {
	var db database.Database

	fname := database.ROOT + "gisaid2020.tsv.gz"
	fmt.Printf("Parsing %s\n", fname)
	db.Parse(fname)

	fmt.Printf("Saving %s\n", database.GOB_NAME)
	db.Save(database.GOB_NAME)

	fmt.Printf("Done\n")
}
