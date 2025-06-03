package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"strings"
)

type CodonFreq struct {
	Codon        string
	CountPer1000 float64
	Aa           byte
	RSCU         float64 // Relative Synonymous Codon Usage

	// The "relative adaptiveness" is like a normalized version of the RSCU.
	// The most preferred codon will end up as 1.0
	RelAd            float64
}

type CodonFreqTable map[string]CodonFreq

// What's the highest RSCU
func (ct CodonFreqTable) Optimum(aa byte) float64 {
	syns := genomes.ReverseCodonTable[aa]
	var best float64
	for _, syn := range syns {
		best = max(best, ct[syn].RSCU)
	}
	return best
}

// Update Aa, RSCU and RelAd for all the CodonFreqs in ct.
func (ct CodonFreqTable) UpdateRSCUs() {
	for k, v := range ct {
		aa := genomes.CodonTable[v.Codon]
		syns := genomes.ReverseCodonTable[aa]

		total := 0.0
		for _, syn := range syns {
			total += ct[syn].CountPer1000
		}

		// What count per 1000 would we expect if all codons synonymous for
		// this AA were used equally?
		expected := total / float64(len(syns))

		// Then the RSCU is just what we do see over that expectation
		rscu := v.CountPer1000 / expected
		w := rscu / ct.Optimum(aa)
		ct[k] = CodonFreq{v.Codon, v.CountPer1000, aa, rscu, w}
	}
}

// Parse a table pasted out of e.g.
// https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=tisspec
func Parse(fname string) CodonFreqTable {
	ret := make(CodonFreqTable)
	utils.Lines(fname, func(line string, err error) bool {
		fields := strings.Split(line, "\t")
		for i := 0; i < len(fields); i += 3 {
			codon := strings.Trim(fields[i], " ")
			if codon == "" {
				continue
			}
			count := utils.Atof(strings.Trim(fields[i+1], " "))
			ret[codon] = CodonFreq{codon, count, 0, 0, 0}
		}
		return true
	})
	ret.UpdateRSCUs()
	return ret
}

type CAICodon struct {
	genomes.Codon
	CAI	float64		// How "preferred" this codon is compared to some reference
}

type CAITranslation []CAICodon

func (c CAITranslation) Mean() float64 {
	var total float64
	for _, codon := range c {
		total += codon.CAI
	}
	return total/float64(len(c))
}

func MakeCAITranslation(g *genomes.Genomes,
	which int, ref CodonFreqTable) CAITranslation {
	trans := genomes.Translate(g, which)
	ret := make(CAITranslation, len(trans))
	for i, codon := range trans {
		cf := ref[codon.Nts]
		best := ref.Optimum(codon.Aa)	// TODO cache these if you care
		cai := cf.RSCU / best
		ret[i] = CAICodon{codon, cai}
	}
	return ret
}

func main() {
	ref := Parse("./human-lung")
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)

	g = g.Filter(0)
	g.Truncate(21562, 25384)

	caiTrans := MakeCAITranslation(g, 0, ref)
	fmt.Println(caiTrans.Mean())
}
