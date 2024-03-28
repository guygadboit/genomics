package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"strings"
)

var orfNames []string

type Alleles map[byte][]int

func joinInts(ints []int, sep string) string {
	s := make([]string, len(ints))
	for i, v := range ints {
		s[i] = fmt.Sprintf("%d", v)
	}
	return strings.Join(s, sep)
}

// Look for alleles where n viruses share one thing, and everybody else has
// the same other thing.
func (a Alleles) checkNearlyUnique(codon genomes.Codon,
	g *genomes.Genomes, n int) {
	if len(a) != 2 {
		return
	}

	var us, them byte
	common := make([]int, n)

	for k, v := range a {
		if k == '-' {
			continue
		}
		if len(v) == n {
			us = k
			copy(common, v)
		} else {
			them = k
		}
	}

	if us != 0 && them != 0 {
		orfs := g.Orfs
		_, orf, pos := codon.OrfRelative(orfs)

		fmt.Printf("%s:%d: %s got %c, everyone else has %c\n",
			orfNames[orf], pos/3+1, joinInts(common, ","), us, them)
	}
}

func main() {
	var numSharers int

	flag.IntVar(&numSharers, "n", 1,
		"Number of genomes sharing same unusual thing")
	flag.Parse()

	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	g.RemoveGaps()

	translations := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		translations[i] = genomes.Translate(g, i)
	}

	for i := 0; i < len(translations[0]); i++ {
		alleles := make(Alleles)
		for j := 0; j < g.NumGenomes(); j++ {
			if i >= len(translations[j]) {
				continue
			}
			codon := translations[j][i]
			aa := codon.Aa
			_, there := alleles[aa]
			if !there {
				alleles[aa] = make([]int, 0)
			}
			alleles[aa] = append(alleles[aa], j)
		}
		alleles.checkNearlyUnique(translations[0][i], g, numSharers)
	}
}

func init() {
	orfNames = []string{
		"ORF1ab",
		"ORF1ab",
		"S",
		"ORF3a",
		"E",
		"M",
		"ORF6",
		"ORF7a",
		"ORF7b",
		"ORF8",
		"N",
		"ORF10",
	}
}
