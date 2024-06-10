package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
)

func main() {
	var db database.Database

	fname := database.ROOT + "gisaid2020.tsv.gz"
	fmt.Printf("Parsing %s\n", fname)
	db.Parse(fname)

	ref := genomes.LoadGenomes("../../fasta/WH1.fasta",
		"../../fasta/WH1.orfs", false)
	db.DetermineSilence(ref)

	fmt.Printf("Saving %s\n", database.GOB_NAME)
	db.Save(database.GOB_NAME)

	fmt.Printf("Done\n")
}
