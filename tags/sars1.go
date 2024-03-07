package main

import (
	"fmt"
	"genomics/genomes"
	"slices"
)

type Similarity struct {
	index int
	SARS1 float64
	SARS2 float64
}

// Print out how close they all are to SARS1 and SARS2, sorted by SARS1
// proximity. Also save the ones that are closer to SARS2
func SARS1Similarity(g *genomes.Genomes) {
	similarities := make([]Similarity, g.NumGenomes())
	SC2Closer := make([]int, 0)

	Tor2 := 79
	SC2 := 0

	for i := 0; i < g.NumGenomes(); i++ {
		ss1 := g.SequenceSimilarity(Tor2, i)
		ss2 := g.SequenceSimilarity(SC2, i)
		similarities[i] = Similarity{i, ss1, ss2}
		if ss2 > ss1 {
			SC2Closer = append(SC2Closer, i)
		}

	}
	slices.SortFunc(similarities, func(a, b Similarity) int {
		if a.SARS1 < b.SARS1 {
			return -1
		} else {
			return 1
		}
	})

	fmt.Printf("%d are closer to sc2\n", len(SC2Closer))

	for _, s := range similarities {
		closer := s.SARS2 > s.SARS1
		fmt.Printf("%d: %f %f %t\n", s.index, s.SARS1*100, s.SARS2*100, closer)
	}

	g.SaveSelected("SARS2-relatives.fasta", SC2Closer...)
}
