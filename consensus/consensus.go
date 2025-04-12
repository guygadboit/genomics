package main

import (
	"fmt"
	"flag"
	"genomics/genomes"
	"genomics/utils"
)

func FindMostSimilar(g *genomes.Genomes, which int, n int) []int {
	type Result struct {
		similarity float64
		index      int
	}
	results := make([]Result, 0, g.NumGenomes()-1)
	for i := 0; i < g.NumGenomes(); i++ {
		if i == which {
			continue
		}
		similarity := g.SequenceSimilarity(which, i)
		r := Result{similarity, i}
		results = append(results, r)
	}

	utils.SortByKey(results, false, func(r Result) float64 {
		return r.similarity
	})

	n = utils.Min(n, len(results))
	ret := make([]int, n)
	for i := 0; i < n; i++ {
		ret[i] = results[i].index
	}
	return ret
}

func Consensus(g *genomes.Genomes, which []int) *genomes.Genomes {
	ret := g.Filter(0, 1)
	ret.DeepCopy(1)

	for i := 0; i < g.Length(); i++ {
		nts := make(map[byte]int)
		for _, w := range which {
			nts[g.Nts[w][i]]++
		}
		var (
			best byte
			most int
		)
		for k, v := range nts {
			if v > most {
				most = v
				best = k
			}
		}
		ret.Nts[1][i] = best
	}
	return ret
}

func main() {
	var fasta string

	flag.StringVar(&fasta, "g", "../fasta/nofcs.fasta", "Alignment to use")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, "../fasta/WH1.orfs", false)

	for i := 0; i < g.NumGenomes(); i++ {
		similar := FindMostSimilar(g, i, 10)
		c := Consensus(g, similar)
		ss := c.SequenceSimilarity(0, 1)
		fmt.Println(ss, i, g.Names[i], similar)
	}
}
