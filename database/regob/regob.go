package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"strings"
)

func main() {
	var (
		fastaName string
		orfs      string
		fname     string
	)

	flag.StringVar(&fastaName, "fasta", "", "Fasta name")
	flag.StringVar(&orfs, "orfs", "../../fasta/WH1.orfs", "ORFs file")
	flag.StringVar(&fname,
		"tsv", database.ROOT+"gisaid2020.tsv.gz", "Include genomes from TSV")
	flag.Parse()

	if fname == "none" {
		fname = ""
	}

	var db database.Database

	if fname != "" {
		fmt.Printf("Parsing %s\n", fname)
		db.Parse(fname)

		ref := genomes.LoadGenomes("../../fasta/WH1.fasta",
			"../../fasta/WH1.orfs", false)
		db.DetermineSilence(ref)
	}

	if fastaName != "" {
		fmt.Printf("Adding in SARS2 relatives from FASTA file\n")
		g := genomes.LoadGenomes(fastaName, orfs, false)

		getHost := func(i int) string {
			if strings.Contains(g.Names[i], "Pangolin") {
				return "Pangolin"
			} else {
				return "Bat"
			}
		}

		database.AddFromGenomes(&db, g, getHost)
	}

	db.BuildMutationIndices()
	db.BuildAccessionIndex()
	db.AddSRAs("read_info2.txt.gz")

	fmt.Printf("Saving %s\n", database.GOB_NAME)
	db.Save(database.GOB_NAME)

	fmt.Printf("Done\n")
}
