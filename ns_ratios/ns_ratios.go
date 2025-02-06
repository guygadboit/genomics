package main

import (
	"fmt"
	"genomics/comparison"
	"genomics/genomes"
)

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives-short-names.fasta",
		"../fasta/WH1.orfs", false)
	g.Truncate(21562, 25384) // spike only

    /*
    // Interestingly none of these have nearly as high a S/N
	g := genomes.LoadGenomes("../fasta/SARS1-relatives.fasta",
		"../fasta/SARS1.orfs", false)
	g.Truncate(21492, 25259) // spike only (SARS1)
    */

	for i := 0; i < g.NumGenomes(); i++ {
		for j := i + 1; j < g.NumGenomes(); j++ {
			c := comparison.Compare(g, i, j)
			S, NS, _ := c.SilentCount()
			ratio := float64(S) / float64(NS)

			fmt.Printf("%.2f %d %d %s vs %s\n", ratio,
				NS, S, g.Names[i], g.Names[j])
		}
	}
}
