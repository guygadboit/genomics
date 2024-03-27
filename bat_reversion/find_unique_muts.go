package main

import (
	"fmt"
	"genomics/genomes"
)

func checkUnique(alleles map[byte][]int,
	codon genomes.Codon, orfs []genomes.Orf) {
	names := []string{
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

	/*
	For SARS1
	names := []string{
		"ORF1ab",
		"ORF1ab",
		"S",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
		"Unknown",
	}
	*/

	for k, v := range alleles {
		if len(v) == 1 && len(alleles) == 2 {
			var other byte
			for k, v := range alleles {
				if len(v) != 1 {
					other = k
					break
				}
			}

			// Don't show indels
			if k == '-' || other == '-' {
				continue
			}

			_, orf, pos := codon.OrfRelative(orfs)
			fmt.Printf("Only in %d: %s: %c%d%c\n",
				v[0], names[orf], k, pos/3+1, other)
		}
	}
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)

	/*
	g := genomes.LoadGenomes("../fasta/SARS1-only-relatives.fasta",
		"../fasta/SARS1.orfs", false)
	*/

	translations := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		translations[i] = genomes.Translate(g, i)
	}

	for i := 0; i < len(translations[0]); i++ {
		alleles := make(map[byte][]int)
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
		checkUnique(alleles, translations[0][i], g.Orfs)
	}
}
