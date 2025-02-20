package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"os"
)

const ROOT = "/fs/f/genomes/viruses/SARS2/"

func main() {
	var reconstruct bool
	var reference string
	var orfs string

	flag.BoolVar(&reconstruct, "reconstruct", false, "Reconstruct fasta files")
	flag.StringVar(&reference, "ref", ROOT+"WH1.fasta", "Reference genome")
	flag.StringVar(&orfs, "orfs", ROOT+"WH1.orfs", "Reference genome ORFs")
	flag.Parse()

	db := database.NewDatabase()

	// EPI_ISL_861438 is an example with some insertions and deletions you can
	// test on.

	var g *genomes.Genomes
	if reconstruct {
		g = genomes.LoadGenomes(reference, orfs, false)
	}

	for _, accNum := range flag.Args() {
		ids := db.GetByAccession(accNum)
		if len(ids) == 0 {
			fmt.Fprintf(os.Stderr, "Can't find %s\n", accNum)
			continue
		}
		id := ids[0]
		r := db.Records[id]

		fmt.Println(r.Summary())

		if reconstruct {
			rg, err := db.Reconstruct(id, g, accNum)
			if err != nil {
				fmt.Fprintln(os.Stderr, err)
				continue
			}
			name := fmt.Sprintf("%s.fasta", r.GisaidAccession)
			rg.SaveMulti(name)
			fmt.Printf("Wrote %s\n", name)
		}
	}
}
