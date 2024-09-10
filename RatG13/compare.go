package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"strings"
    "math"
)

func restrict(g *genomes.Genomes, start, end int) {
	for i := 0; i < g.NumGenomes(); i++ {
		g.Nts[i] = g.Nts[i][start:end]
	}
}

func Compare4991() {
	g := genomes.LoadGenomes("WH1-RaTG13.fasta", "../fasta/WH1.orfs", false)
	f4991 := genomes.LoadGenomes("4991.fasta", "", false)
	n := f4991.Length()
	// fmt.Printf("4991 is %d nts long\n", n)

	for i := 0; i < g.Length()-n; i++ {
		g2 := g.Filter(0, 1)
		restrict(g2, i, i+n)
		fmt.Printf("%d %f\n", i, g2.SequenceSimilarity(0, 1))
	}
}

func ComparePair(g *genomes.Genomes, a, b int) (int, int) {
	boundary := 681 // between S1 and S2
	before, after := 0, 0

	for i := 0; i < boundary; i++ {
		if g.Nts[a][i] != g.Nts[b][i] {
			before++
		}
	}
	for i := boundary; i < g.Length(); i++ {
		if g.Nts[a][i] != g.Nts[b][i] {
			after++
		}
	}
	return before, after
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

func CompareS1S2() {
	g := genomes.LoadGenomes("./spikes.fasta", "", false)
	// g.Names = LoadShortNames()
	for a := 0; a < g.NumGenomes(); a++ {
		for b := 0; b < a; b++ {
			before, after := ComparePair(g, a, b)
            ratio := float64(before)/float64(after)
            if math.IsNaN(ratio) {
                continue
            }

			fmt.Printf("%f %s vs %s. S1: %d changes S2: %d changes\n",
				float64(before)/float64(after),
				g.Names[a], g.Names[b],
				before, after)
		}
	}
}

func main() {
	CompareS1S2()
}
