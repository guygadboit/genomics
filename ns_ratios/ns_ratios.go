package main

import (
	"flag"
	"fmt"
	"genomics/comparison"
	"genomics/genomes"
	"log"
)

func main() {
	var (
		fasta, orfs string
		spikeOnly   bool
	)

	flag.StringVar(&fasta, "fasta",
		"../fasta/SARS2-relatives-short-names.fasta", "Fasta file to use")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "ORFs file to use")
	flag.BoolVar(&spikeOnly, "spike", false, "S only")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	if spikeOnly {
		S, err := g.Orfs.Find("S")
		if err != nil {
			log.Fatal("Can't find S")
		}
		g.Truncate(S.Start, S.End)
	}

	for i := 0; i < g.NumGenomes(); i++ {
		for j := i + 1; j < g.NumGenomes(); j++ {
			c := comparison.Compare(g, i, j)
			S, NS, _ := c.SilentCount()
			ratio := float64(S) / float64(NS)

			fmt.Printf("%.2f %d %d %s(%d) vs %s(%d)\n", ratio,
				NS, S, g.Names[i], i, g.Names[j], j)
		}
	}
}
