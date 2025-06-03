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
	RSCU         float64 // Relative Synonymous Codon Usage
}

type CodonFreqTable map[string]CodonFreq

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
		ct[k] = CodonFreq{v.Codon, v.CountPer1000, rscu}
	}
}

// Parse a table pasted out of e.g.
// https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=tisspec
func Parse(fname string) CodonFreqTable {
	ret := make(CodonFreqTable)
	utils.Lines(fname, func(line string, err error) bool {
		fields := strings.Split(line, "\t")
		fmt.Println(fields)
		for i := 0; i < len(fields); i += 3 {
			codon := strings.Trim(fields[i], " ")
			if codon == "" {
				continue
			}
			count := utils.Atof(strings.Trim(fields[i+1], " "))
			ret[codon] = CodonFreq{codon, count, 0}
		}
		return true
	})
	ret.UpdateRSCUs()
	return ret
}

func main() {
	cft := Parse("./lung")
	for _, v := range cft {
		fmt.Println(v)
	}
}
