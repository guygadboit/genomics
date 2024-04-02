package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"strings"
)

type Alleles map[byte][]int

func joinInts(ints []int, sep string) string {
	s := make([]string, len(ints))
	for i, v := range ints {
		s[i] = fmt.Sprintf("%d", v)
	}
	return strings.Join(s, sep)
}

// Counts how many "quirks" (lone differences from the crowd) in each genome.
// Maps a genome index to a count
type QuirkMap map[int]int

func (q QuirkMap) Combine(other QuirkMap) {
	for k, v := range other {
		q[k] += v
	}
}

// Look for alleles where n viruses share one thing, and everybody else has
// the same other thing.
func (a Alleles) checkNearlyUnique(codon genomes.Codon,
	g *genomes.Genomes, n int) QuirkMap {
	ret := make(QuirkMap)

	if len(a) != 2 {
		return ret
	}

	var us, them byte
	common := make([]int, n)

	for k, v := range a {
		if k == '-' {
			continue
		}
		if len(v) <= n {
			us = k
			copy(common, v)
		} else {
			them = k
		}
	}

	if us != 0 && them != 0 {
		orfs := g.Orfs
		orf, pos, _ := orfs.GetOrfRelative(codon.Pos)

		fmt.Printf("%s:%d: %s got %c, everyone else has %c\n",
		orfs[orf].Name, pos/3+1, joinInts(common, ","), us, them)

		for _, v := range common {
			ret[v] += 1
		}
	}
	return ret
}

// Look for alleles unique to which, but not caring whether the others all have
// the same thing there or are variously different.
func (a Alleles) checkUnique2(codon genomes.Codon,
	g *genomes.Genomes, which int) {
	orfs := g.Orfs

	var us byte
	for k, v := range a {
		if len(v) == 1 && v[0] == which {
			us = k
			break
		}
	}
	if us == 0 {
		return
	}

	orf, pos, _ := orfs.GetOrfRelative(codon.Pos)
	fmt.Printf("%d has %s:%d%c, everyone else something else\n",
		which, orfs[orf].Name, pos/3+1, us)
}

func graphData(qm QuirkMap, g *genomes.Genomes) {
	for k, v := range qm {
		ss := g.SequenceSimilarity(0, k)
		fmt.Println(ss, v)
	}
}

func main() {
	var numSharers int
	var unique int

	flag.IntVar(&numSharers, "n", 1,
		"Number of genomes sharing same unusual thing")
	flag.IntVar(&unique, "u", -1,
		"Just check for unique whatever the others have")
	flag.Parse()

	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	g.RemoveGaps()

	// quirks := make(QuirkMap)

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
		// q := alleles.checkNearlyUnique(translations[0][i], g, numSharers)
		// quirks.Combine(q)

		if unique != -1 {
			alleles.checkUnique2(translations[0][i], g, unique)
		} else {
			alleles.checkNearlyUnique(translations[0][i], g, numSharers)
		}
	}
	// graphData(quirks, g)
}
