package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/degeneracy"
)

type Row struct {
	FourA	float64		// number of codons with 4-fold degeneracy and A in the 3rd pos
	FourC	float64		// etc...
	FourG	float64
	FourT	float64

	TwoA	float64
	TwoC	float64
	TwoG	float64
	TwoT	float64
}

func countClasses(t degeneracy.Translation) Row {
	var ret Row
	var totalFours, totalTwos float64

	for _, c := range t {
		nt := c.Nts[2]
		fold := c.Fold
		if fold == 3 {	// Treat 3-fold as if it were 2
			fold = 2
		}
		switch fold {
		case 4:
			switch nt {
			case 'A':
				ret.FourA++
			case 'G':
				ret.FourG++
			case 'C':
				ret.FourC++
			case 'T':
				ret.FourT++
			}
			totalFours++
		case 2:
			switch nt {
			case 'A':
				ret.TwoA++
			case 'G':
				ret.TwoG++
			case 'C':
				ret.TwoC++
			case 'T':
				ret.TwoT++
			}
			totalTwos++
		}
	}

	ret.FourA /= totalFours
	ret.FourG /= totalFours
	ret.FourC /= totalFours
	ret.FourT /= totalFours

	ret.TwoA /= totalTwos
	ret.TwoG /= totalTwos
	ret.TwoC /= totalTwos
	ret.TwoT /= totalTwos

	return ret
}

func main() {
	g := genomes.LoadGenomes("../../fasta/SARS2-relatives.fasta",
		"../../fasta/WH1.orfs", false)

	for i := 0; i < g.NumGenomes(); i++ {
		t := degeneracy.Translate(g, i)
		classes := countClasses(t)
		fmt.Println(classes)
	}
}
