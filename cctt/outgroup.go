package main

import (
	"fmt"
	"genomics/genomes"
	"slices"
	"strings"
)

type Outgroup struct {
	g *genomes.Genomes
}

func NewOutgroup() *Outgroup {
	g := genomes.LoadGenomes(
		"../fasta/SARS2-relatives-short-names.fasta",
		"../fasta/WH1.orfs", false)
	// g = g.Filter(7, 8, 6, 5)
	return &Outgroup{g}
}

type Alt struct {
	nt    byte
	count int
}

// Sorted by most frequent first
type Alts []Alt

func (al Alts) ToString() string {
	s := make([]string, len(al))
	for i, a := range al {
		s[i] = fmt.Sprintf("%cx%d", a.nt, a.count)
	}
	return strings.Join(s, ",")
}

func (o *Outgroup) Get(pos int) Alts {
	counts := make(map[byte]int)
	for i := 1; i < o.g.NumGenomes(); i++ {
		counts[o.g.Nts[i][pos]]++
	}
	ret := make(Alts, 0, len(counts))
	for k, v := range counts {
		ret = append(ret, Alt{k, v})
	}
	slices.SortFunc(ret, func(a, b Alt) int {
		if a.count > b.count {
			return -1
		}
		if a.count < b.count {
			return 1
		}
		return 0
	})
	return ret
}

func (o *Outgroup) WhoHas(pos int, nt byte) []int {
	g := o.g
	ret := make([]int, 0)
	for i := 1; i < g.NumGenomes(); i++ {
		if g.Nts[i][pos] == nt {
			ret = append(ret, i)
		}
	}
	return ret
}

func (o *Outgroup) DisplaySorted(counts map[int]int) {
	type result struct {
		index, count int
	}
	results := make([]result, 0, len(counts))
	for k, v := range counts {
		results = append(results, result{k, v})
	}
	slices.SortFunc(results, func(a, b result) int {
		if a.count > b.count {
			return -1
		}
		if a.count < b.count {
			return 1
		}
		return 0
	})
	for _, result := range results {
		ss := fmt.Sprintf("%.2f%%", o.g.SequenceSimilarity(0, result.index)*100)
		fmt.Println(result.index, o.g.Names[result.index], result.count, ss)
	}
}
