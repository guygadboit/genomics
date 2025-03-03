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
	var (
		reconstruct bool
		msa         bool
		reference   string
		orfs        string
	)

	flag.BoolVar(&reconstruct, "reconstruct", false, "Reconstruct fasta files")
	flag.BoolVar(&msa, "msa", false, "Output reconstruction of all inputs"+
		"as a single multi-sequence alignment")
	flag.StringVar(&reference, "ref", ROOT+"WH1.fasta", "Reference genome")
	flag.StringVar(&orfs, "orfs", ROOT+"WH1.orfs", "Reference genome ORFs")

	flag.Parse()

	db := database.NewDatabase()

	// EPI_ISL_861438 is an example with some insertions and deletions you can
	// test on.

	var g, output *genomes.Genomes
	if reconstruct {
		g = genomes.LoadGenomes(reference, orfs, false)
		if msa {
			output = g.Filter(0)
		}
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

			if !msa {
				name := fmt.Sprintf("%s.fasta", r.GisaidAccession)
				rg.SaveMulti(name)
				fmt.Printf("Wrote %s\n", name)
				continue
			}

			if len(rg.Nts[1]) != len(g.Nts[0]) {
				continue
			}

			output.Nts = append(output.Nts, rg.Nts[1])
			output.Names = append(output.Names, rg.Names[1])
		}
	}

	if output != nil {
		name := "GISAID-genomes.fasta"
		output.SaveMulti(name)
		fmt.Printf("Wrote %s\n", name)
	}
}
