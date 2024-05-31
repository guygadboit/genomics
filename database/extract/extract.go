package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"os"
	"time"
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

	// Find something with some insertions and deletions
	/*
		for _, r := range db.Records {
			if len(r.Insertions) > 0 &&
				len(r.Deletions) > 0 &&
				len(r.NucleotideChanges) > 0 {
				fmt.Println(r.GisaidAccession, len(r.Insertions), len(r.Deletions))
			}
		}
		return
	*/
	// EPI_ISL_861438

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

		fmt.Printf("%s %s %s %s %s: %s %s; %d %d\n",
			r.GisaidAccession,
			r.CollectionDate.Format(time.DateOnly),
			r.Country,
			r.Region,
			r.City,
			r.NucleotideChanges.ToString(),
			r.AAChanges.ToString(),
			len(r.Insertions),
			len(r.Deletions))

		if reconstruct {
			rg, err := db.Reconstruct(id, g)
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
