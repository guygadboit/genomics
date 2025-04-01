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

type Allele struct {
	nt    byte
	count int
}

// Sorted by most frequent first
type Alleles []Allele

func (al Alleles) ToString() string {
	s := make([]string, len(al))
	for i, a := range al {
		s[i] = fmt.Sprintf("%cx%d", a.nt, a.count)
	}
	return strings.Join(s, ",")
}

func (o *Outgroup) Get(pos int) Alleles {
	counts := make(map[byte]int)
	for i := 1; i < o.g.NumGenomes(); i++ {
		counts[o.g.Nts[i][pos]]++
	}
	ret := make(Alleles, 0, len(counts))
	for k, v := range counts {
		ret = append(ret, Allele{k, v})
	}
	slices.SortFunc(ret, func(a, b Allele) int {
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
