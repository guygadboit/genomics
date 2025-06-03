package main

import (
	"fmt"
	"strings"
	"genomics/genomes"
	"genomics/utils"
)

type CodonFreq struct {
	Codon        string
	CountPer1000 float64
	RSCU         float64 // Relative Synonymous Codon Usage
}

type CodonFreqTable map[string]CodonFreq

func (c *CodonFreq) UpdateRSCU() {
	aa := genomes.CodonTable[c.Codon]
	numSyn := len(genomes.ReverseCodonTable[aa])


	// Wrong. That's how often this codon appears, not how often out of
	// possible synonyms for the same one. But you do have that information.
	count := c.CountPer1000 / 1000
	c.RSCU = count / float64(numSyn)
}

func (c *CodonFreq) Init(codon string, count1000 float64) {
	c.Codon = codon
	c.CountPer1000 = count1000
	c.UpdateRSCU()
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
			var cf CodonFreq
			cf.Init(codon, count)
			ret[cf.Codon] = cf
		}
		return true
	})
	return ret
}

func main() {
	cft := Parse("./lung")
	fmt.Println(cft)
}
