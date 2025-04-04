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
	return &Outgroup{genomes.LoadGenomes(
		"../fasta/SARS2-relatives-short-names.fasta",
		"../fasta/WH1.orfs", false)}
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
