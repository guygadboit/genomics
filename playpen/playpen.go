package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
)

type CodonMap map[byte][]genomes.Codon

func (c CodonMap) Add(codon genomes.Codon) {
	key := codon.Aa
	if _, there := c[key]; !there {
		c[key] = make([]genomes.Codon, 0)
	}
	c[key] = append(c[key], codon)

}

func (c CodonMap) Intersection(other CodonMap) CodonMap {
	ret := make(CodonMap)
	for k, v := range c {
		if _, there := other[k]; there {
			ret[k] = append(v, other[k]...)
		}
	}
	return ret
}

func CountQN(g *genomes.Genomes) {
	count := func(seq []byte) int {
		trans := genomes.TranslateAlignedShort(seq)
		ret := 0
		for _, c := range trans {
			if c == 'Q' || c == 'N' {
				ret++
			}
		}
		fmt.Println(ret, string(trans))
		return ret
	}

	l := 60 * 3
	numSeqs := 0
	for i := 0; i < g.Length()-l+1; i++ {
		seq := g.Nts[0][i : i+l]
		count(seq)
		count(utils.ReverseComplement(seq))
		numSeqs += 2
	}
	fmt.Println("num seqs", numSeqs)
}

func ConvergentPangolin(g *genomes.Genomes) {
	gx := utils.ToSet([]int{35, 36, 37, 38, 39, 40, 41})
	gd := utils.ToSet([]int{30, 31, 32})

	translations := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		translations[i] = genomes.Translate(g, i)
	}

	numCodons := len(translations[0])

	for i := 0; i < numCodons; i++ {
		gxCodons := make(CodonMap)
		gdCodons := make(CodonMap)
		otherCodons := make(CodonMap)

		for j, trans := range translations {
			codon := trans[i]
			if codon.Aa == '-' {
				continue
			}
			if gx[j] {
				gxCodons.Add(codon)
			} else if gd[j] {
				gdCodons.Add(codon)
			} else {
				otherCodons.Add(codon)
			}
		}

		// ix := gxCodons
		ix := gdCodons.Intersection(gxCodons)

		for k, v := range ix {
			if _, there := otherCodons[k]; there {
				continue
			}
			orfs := g.Orfs
			for _, c := range v {
				orfI, oPos, _ := orfs.GetOrfRelative(c.Pos)
				name := orfs[orfI].Name
				fmt.Printf("%s:%d%c\n", name, oPos/3+1, k)
			}
		}
	}
}

func main() {
	if false {
		g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
			"../fasta/WH1.orfs", false)
		CountQN(g)
	}

	if false {
		config := stats.BlastDefaultConfig()
		br := stats.Blast(config, "viruses/SARS2", []byte("TGGTCGC"),
			10, 10, stats.VERBOSE)
		fmt.Println(br)
	}

	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	muts := mutations.PossibleSilentMuts(g, 0)

	for _, mut := range muts {
		if mut.Pos >= 22561 && mut.Pos < 23600 {
			fmt.Println(mut)
		}
	}
	// ConvergentPangolin(g)
}
