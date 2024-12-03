package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/genomes"
)

func main() {
	var fromFasta bool

	flag.BoolVar(&fromFasta, "fasta", false, "Include genomes from fasta")
	flag.Parse()

	var db database.Database

	fname := database.ROOT + "gisaid2020.tsv.gz"
	fmt.Printf("Parsing %s\n", fname)
	db.Parse(fname)

	ref := genomes.LoadGenomes("../../fasta/WH1.fasta",
		"../../fasta/WH1.orfs", false)
	db.DetermineSilence(ref)

	if fromFasta {
		fmt.Printf("Adding in SARS2 relatives from FASTA file\n")
		g := genomes.LoadGenomes("../../fasta/SARS2-relatives.fasta",
			"../../fasta/WH1.orfs", false)
		database.AddFromGenomes(&db, g, nil)
	}

	db.BuildMutationIndices()
	db.BuildAccessionIndex()

	fmt.Printf("Saving %s\n", database.GOB_NAME)
	db.Save(database.GOB_NAME)

	fmt.Printf("Done\n")
}
