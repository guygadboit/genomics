package main

import (
	"flag"
	"fmt"
	"genomics/comparison"
	"genomics/genomes"
	"genomics/mutations"
	"log"
)

func ShowPossible(g *genomes.Genomes) {
	fmt.Println("Possible silent muts")
	ct := 0
	for i := 0; i < g.NumGenomes(); i++ {
		muts := mutations.PossibleSilentMuts(g, i)
		for _, m := range muts {
			if m.From == 'C' && m.To == 'T' {
				ct++
			}
		}
		fmt.Printf("%d of which %d are CT %s\n", len(muts), ct, g.Names[i])
	}
}

func main() {
	var (
		fasta, orfs string
		spikeOnly   bool
		possible    bool
	)

	flag.StringVar(&fasta, "fasta",
		"../fasta/SARS2-relatives-short-names.fasta", "Fasta file to use")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "ORFs file to use")
	flag.BoolVar(&spikeOnly, "spike", false, "S only")
	flag.BoolVar(&possible, "poss", false, "Show possible")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)
	g.RemoveGaps()

	if spikeOnly {
		S, err := g.Orfs.Find("S")
		if err != nil {
			log.Fatal("Can't find S")
		}
		g.Truncate(S.Start, S.End)
	}

	if possible {
		ShowPossible(g)
		return
	}

	for i := 0; i < g.NumGenomes(); i++ {
		for j := i + 1; j < g.NumGenomes(); j++ {
			c := comparison.Compare(g, i, j)
			S, NS, _ := c.SilentCount()
			ratio := float64(S) / float64(NS)

			ss := g.SequenceSimilarity(i, j, false)
			fmt.Printf("%.2f %d %d %s(%d) vs %s(%d) %.4f SS\n", ratio,
				NS, S, g.Names[i], i, g.Names[j], j, ss)
		}
	}
}
