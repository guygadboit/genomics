package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"math"
	"strings"
)

type CodonFreq struct {
	Codon        string
	CountPer1000 float64
	Aa           byte
	RSCU         float64 // Relative Synonymous Codon Usage

	// The "relative adaptiveness" is like a normalized version of the RSCU.
	// The most preferred codon will end up as 1.0
	RelAd float64
}

type CodonFreqTable struct {
	Frequencies map[string]CodonFreq
	Optima      map[byte]float64 // The best RSCU for each AA
}

// What's the highest RSCU for a particular AA?
func (ct CodonFreqTable) CalcOptimum(aa byte) float64 {
	syns := genomes.ReverseCodonTable[aa]
	var best float64
	for _, syn := range syns {
		best = max(best, ct.Frequencies[syn].RSCU)
	}
	return best
}

func (ct *CodonFreqTable) UpdateOptima() {
	ct.Optima = make(map[byte]float64)
	for k, _ := range genomes.ReverseCodonTable {
		ct.Optima[k] = ct.CalcOptimum(k)
	}
}

// Update Optima and Aa, RSCU and RelAd for all the CodonFreqs in
// ct.Frequencies
func (ct *CodonFreqTable) UpdateComputedValues() {
	for k, v := range ct.Frequencies {
		aa := genomes.CodonTable[v.Codon]
		syns := genomes.ReverseCodonTable[aa]

		total := 0.0
		for _, syn := range syns {
			total += ct.Frequencies[syn].CountPer1000
		}

		// What count per 1000 would we expect if all codons synonymous for
		// this AA were used equally?
		expected := total / float64(len(syns))

		// Then the RSCU is just what we do see over that expectation
		rscu := v.CountPer1000 / expected
		ct.Frequencies[k] = CodonFreq{v.Codon, v.CountPer1000, aa, rscu, 0}
	}

	// Once we know all the RSCUs we can work out the relative adaptations
	ct.UpdateOptima()

	// And now we can work out the relative adaptations
	for k, v := range ct.Frequencies {
		newCf := v
		newCf.RelAd = v.RSCU / ct.Optima[v.Aa]
		ct.Frequencies[k] = newCf
	}
}

// Parse a table pasted out of e.g.
// https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=tisspec
func Parse(fname string) *CodonFreqTable {
	ret := new(CodonFreqTable)
	ret.Frequencies = make(map[string]CodonFreq)
	utils.Lines(fname, func(line string, err error) bool {
		fields := strings.Split(line, "\t")
		for i := 0; i < len(fields); i += 3 {
			codon := strings.Trim(fields[i], " ")
			if codon == "" {
				continue
			}
			count := utils.Atof(strings.Trim(fields[i+1], " "))
			ret.Frequencies[codon] = CodonFreq{codon, count, 0, 0, 0}
		}
		return true
	})
	ret.UpdateComputedValues()
	return ret
}

type RelCodon struct {
	genomes.Codon
	RelAd float64 // The relative adaptation of this codon
}

// A translation decorated with relative adaptation values
type RelTranslation []RelCodon

func MakeRelTranslation(g *genomes.Genomes,
	which int, ref *CodonFreqTable) (RelTranslation, float64) {
	trans := genomes.Translate(g, which)
	ret := make(RelTranslation, len(trans))

	prodRSCU, prodBest := 1.0, 1.0
	exceptions := 0

	for i, codon := range trans {
		cf := ref.Frequencies[codon.Nts]
		best := ref.Optima[codon.Aa]
		ret[i] = RelCodon{codon, cf.RelAd}

		prodRSCU *= cf.RSCU
		prodBest *= best

		switch codon.Nts {
		case "ATG":
			fallthrough
		case "TGG":
			exceptions++
		}
	}

	l := float64(len(trans) - exceptions)
	caiObs := math.Pow(prodRSCU, 1/l)
	caiMax := math.Pow(prodBest, 1/l)

	return ret, caiObs / caiMax
}

func main() {
	ref := Parse("./human-lung")
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)

	g = g.Filter(0)
	start := 21562
	end := 25384
	g.Truncate(start, end)
	g.SaveMulti("check.fasta")

	// g := genomes.LoadGenomes("check.fasta", "", false)

	relTrans, cai := MakeRelTranslation(g, 0, ref)
	fmt.Println("CAI", cai)
	fmt.Println(ref.Optima['L'])

	for _, rc := range relTrans {
		fmt.Printf("%s %c %f\n", rc.Nts, rc.Aa, rc.RelAd)
	}
}
