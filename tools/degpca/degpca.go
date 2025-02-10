package main

import (
	"fmt"
	"genomics/degeneracy"
	"genomics/genomes"
	"genomics/stats"
)

/*
Contains 8 items, the proportion of codons with 4-fold degeneracy that have A,
C, G, and T in the 3rd position, followed by the same thing for 2-fold
*/
type Row []float64

func countClasses(t degeneracy.Translation) Row {
	ret := make(Row, 8)
	var totalFours, totalTwos float64

	for _, c := range t {
		nt := c.Nts[2]
		fold := c.Fold
		if fold == 3 { // Treat 3-fold as if it were 2
			fold = 2
		}

		var i int
		if fold == 4 {
			i = 0
			totalFours++
		} else {
			i = 4
			totalTwos++
		}

		var j int
		switch nt {
		case 'A':
			j = 0
		case 'G':
			j = 1
		case 'C':
			j = 2
		case 'T':
			j = 3
		}

		ret[i+j]++
	}

	for i := 0; i < 4; i++ {
		ret[i] /= totalFours
	}
	for i := 4; i < 8; i++ {
		ret[i] /= totalTwos
	}

	return ret
}

func main() {
	g := genomes.LoadGenomes("../../fasta/SARS2-relatives.fasta",
		"../../fasta/WH1.orfs", false)

	data := make([][]float64, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		t := degeneracy.Translate(g, i)
		classes := countClasses(t)
		data = append(data, classes)
	}

	result := stats.PCA(2, data)
	fmt.Println(result)
}
