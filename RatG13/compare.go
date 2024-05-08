package main

import (
	"genomics/genomes"
	"fmt"
)

func restrict(g *genomes.Genomes, start, end int) {
	for i := 0; i < g.NumGenomes(); i++ {
		g.Nts[i] = g.Nts[i][start:end]
	}
}

func main() {
	g := genomes.LoadGenomes("WH1-RaTG13.fasta", "../fasta/WH1.orfs", false)
	f4991 := genomes.LoadGenomes("4991.fasta", "", false)
	n := f4991.Length()
	// fmt.Printf("4991 is %d nts long\n", n)

	for i := 0; i < g.Length() - n; i++ {
		g2 := g.Filter(0, 1)
		restrict(g2, i, i+n)
		fmt.Printf("%d %f\n", i, g2.SequenceSimilarity(0, 1))
	}
}
