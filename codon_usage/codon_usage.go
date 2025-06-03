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

	/*
	The CAI is relative to the CodonFreq for some other organism. For example,
	you might generate a CodonFreqTable for human lung, and then another one
	for some virus, and set the CAIs in the virus one to tell you how well its
	codon usage matched that of the human lung.
	*/
	CAI          float64 // Codon Adaptation Index
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
		ct[k] = CodonFreq{v.Codon, v.CountPer1000, aa, rscu, 0}
	}
}

// What's the highest RSCU we have for this AA?
func (ct CodonFreqTable) Optimum(aa byte) float64 {
	syns := genomes.ReverseCodonTable[aa]
	var best float64
	for _, syn := range syns {
		if ct[syn].RSCU > best {
			best = ct[syn].RSCU
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
			ret[codon] = CodonFreq{codon, count, 0, 0, 0}
		}
		return true
	})
	ret.UpdateRSCUs()
	return ret
}

func FindCAITable(g *genomes.Genomes,
	which int, ref CodonFreqTable) (CodonFreqTable, float64) {
	ret := make(CodonFreqTable)
	trans := genomes.Translate(g, which)
	counts := make(map[string]int)
	for _, codon := range trans {
		counts[codon.Nts]++
	}
	for k, v := range counts {
		ret[k] = CodonFreq{k, float64(v) / 1000, 0, 0, 0}
	}

	ret.UpdateRSCUs()

	total := 0.0
	for k, v := range ret {
		opt := ref.Optimum(v.Aa)

		v.CAI = v.RSCU / opt
		ret[k] = v
		total += v.CAI
	}

	return ret, total/float64(len(ret))
}

func main() {
	ref := Parse("./lung")
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	caiTable, mean := FindCAITable(g, 0, ref)

	for _, v := range caiTable {
		fmt.Printf("%s %f %c %f\n", v.Codon, v.RSCU, v.Aa, v.CAI)
	}
	fmt.Printf("Mean: %f\n", mean)
}
