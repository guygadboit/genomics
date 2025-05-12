package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
	"log"
)

func Match(g *genomes.Genomes, pu *pileup.Pileup, ref int, verbose bool) int {
	var ret int

	for i := 0; i < g.Length(); i++ {

		// Skip the actual FCS
		if i >= 23600 && i < 23612 {
			continue
		}

		rec := pu.Get(i)
		if rec == nil {
			continue
		}
		p2s := rec.Reads[0].Nt
		refNt := g.Nts[ref][i]

		others := 0
		for j := 0; j < g.NumGenomes(); j++ {
			if j == ref {
				continue
			}
			if g.Nts[j][i] == p2s {
				others++
			}
		}

		if others == 0 {
			var same string
			if refNt == p2s {
				same = " SAME"
				ret++
			}
			if verbose {
				fmt.Printf("%d: %c %c%s\n", i+1, refNt, p2s, same)
			}
		}
	}
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)

	pu, err := pileup.Parse("./P2S-pileup")
	if err != nil {
		log.Fatal(err)
	}

	for i := 0; i < g.NumGenomes(); i++ {
		fmt.Println(i, Match(g, pu, i, false))
	}
}
