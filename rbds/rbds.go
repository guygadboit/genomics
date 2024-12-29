package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"slices"
	"strings"
)

var SPIKE_START utils.OneBasedPos = utils.OneBasedPos(21563)

type Comparison struct {
	index      int
	similarity float64
}

func printComparisons(g *genomes.Genomes, c []Comparison, caption string) {
	slices.SortFunc(c, func(a, b Comparison) int {
		if a.similarity < b.similarity {
			return 1
		}
		if a.similarity > b.similarity {
			return -1
		}
		return 0
	})
	for _, c := range c {
		fmt.Printf("%s %s %.4f\n", caption, g.Names[c.index], c.similarity)
	}
}

func Compare(g *genomes.Genomes, start, end utils.OneBasedPos) {
	comparisons := make([]Comparison, 0)
	s := int((SPIKE_START + start*3) - 1)
	e := int(SPIKE_START + end*3)
	g.Truncate(s, e)

	// Protein first
	refTrans := genomes.Translate(g, 0)

	for i := 1; i < g.NumGenomes(); i++ {
		same := 0
		trans := genomes.Translate(g, i)
		for j, codon := range refTrans {
			if codon.Aa == trans[j].Aa {
				same++
			}
		}
		similarity := float64(same) / float64(len(refTrans))
		comparisons = append(comparisons, Comparison{i, similarity})
	}
	printComparisons(g, comparisons, "Protein")

	comparisons = make([]Comparison, 0)
	for i := 1; i < g.NumGenomes(); i++ {
		same := 0
		for j := 0; j < g.Length(); j++ {
			if g.Nts[0][j] == g.Nts[i][j] {
				same++
			}
		}
		similarity := float64(same) / float64(g.Length())
		comparisons = append(comparisons, Comparison{i, similarity})
	}
	printComparisons(g, comparisons, "Nucleotide")
}

func LoadShortNames() []string {
	ret := make([]string, 0)
	utils.Lines("../fasta/short_names.txt", func(line string, err error) bool {
		fields := strings.Split(line, " ")
		ret = append(ret, fields[1])
		return true
	})
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	g.Names = LoadShortNames()

	// As AA offsets into S
	rbd := []utils.OneBasedPos{333, 679}
	rbm := []utils.OneBasedPos{438, 506}

	fmt.Println("RBD")
	Compare(g.Clone(), rbd[0], rbd[1])
	fmt.Println("RBM")
	Compare(g.Clone(), rbm[0], rbm[1])
}
