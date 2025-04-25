package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

func ContigMatches(g *genomes.Genomes, a, b int) map[int]int {
	ret := make(map[int]int)
	for start := 0; start < g.Length(); start++ {
		length := 0
		for i := start; i < g.Length(); i++ {
			if g.Nts[a][i] != g.Nts[b][i] {
				break
			}
			length++
		}
		ret[length]++
	}
	return ret
}

func DisplayMatches(matches map[int]int) {
	type result struct {
		length int
		count int
	}
	results := make([]result, 0, len(matches))
	for k, v := range matches {
		results = append(results, result{k, v})
	}
	utils.SortByKey(results, false, func(r result) int {
		return r.length
	})
	for _, r := range results {
		fmt.Printf("%d: %d\n", r.length, r.count)
	}
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	matches := ContigMatches(g, 7, 8)
	DisplayMatches(matches)
}
