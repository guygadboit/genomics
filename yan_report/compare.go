package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"slices"
	"strings"
)

type Result struct {
	which  int
	muts   float64 // as a fraction
	nsMuts float64
}

/* Returns the sequence similarity and the protein similarity */
func CompareOrf(g *genomes.Genomes,
	a, b int, orf genomes.Orf) Result {
	var muts, nsMuts int
	for i := orf.Start; i < orf.End; i++ {
		if g.Nts[a][i] == g.Nts[b][i] {
			continue
		}
		silent, _, _ := genomes.IsSilent(g, i, 1, a, b)
		muts++

		if !silent {
			nsMuts++
		}
	}
	total := float64(orf.End - orf.Start)
	mutsF, nsMutsF := (total-float64(muts))/total, (total-float64(nsMuts))/total
	return Result{b, mutsF, nsMutsF}
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

func CompareAll(g *genomes.Genomes) {
	shortNames := LoadShortNames()
	for _, orf := range g.Orfs {
		results := make([]Result, 0)
		for i := 1; i < g.NumGenomes(); i++ {
			result := CompareOrf(g, 0, i, orf)
			results = append(results, result)
		}

		slices.SortFunc(results, func(a, b Result) int {
			if a.nsMuts < b.nsMuts {
				return 1
			}
			if a.nsMuts > b.nsMuts {
				return -1
			}
			return 0
		})

		fmt.Println(orf.Name)
		for _, result := range results {
			fmt.Printf("%-16s: %s nt=%.2f, prot=%.2f\n",
				shortNames[result.which], orf.Name, result.muts, result.nsMuts)
		}
		fmt.Println()
	}
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	CompareAll(g)
}
