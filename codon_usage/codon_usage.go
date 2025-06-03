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
}

type CodonFreqTable map[string]CodonFreq

// Update Aa and RSCU for all the CodonFreqs in ct.
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
		ct[k] = CodonFreq{v.Codon, v.CountPer1000, aa, rscu}
	}
}

func (ct CodonFreqTable) Optimum(aa byte) CodonFreq {
	syns := genomes.ReverseCodonTable[aa]
	var best CodonFreq
	for _, syn := range syns {
		if ct[syn].RSCU > best.RSCU {
			best = ct[syn]
		}
	}
	return best
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
			ret[codon] = CodonFreq{codon, count, 0, 0}
		}
		return true
	})
	ret.UpdateRSCUs()
	return ret
}

func FindFreq(g *genomes.Genomes,
	which int, ref CodonFreqTable) CodonFreqTable {
	ret := make(CodonFreqTable)
	trans := genomes.Translate(g, which)
	counts := make(map[string]int)
	for _, codon := range trans {
		counts[codon.Nts]++
	}
	for k, v := range counts {
		ret[k] = CodonFreq{k, float64(v)/1000, 0, 0}
	}
	ret.UpdateRSCUs()
	return ret
}

func main() {
	ref := Parse("./lung")
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	cft := FindFreq(g, 0, ref)

	for _, v := range cft {
		opt := ref.Optimum(v.Aa)
		fmt.Printf("%s %f %c %f\n", v.Codon, v.RSCU, v.Aa, v.RSCU/opt)
	}
}
