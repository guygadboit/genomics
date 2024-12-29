package main

import (
	"fmt"
	"strings"
	"genomics/genomes"
	"genomics/utils"
)

var BASIC map[byte]bool

func init() {
	BASIC = utils.ToSet([]byte("KRH"))
}

type NLS struct {
	codons	[]genomes.Codon		// 7 of these
}

func (n *NLS) IsValid() bool {
	if n.codons[0].Aa != 'P' {
		return false
	}

	// Now somewhere in the next 6 we need a run of 3 out of 4 that are basic
	for i := 0; i < 2; i++ {
		count := 0
		for j := i; j < i+4; j++ {
			if BASIC[n.codons[j].Aa] {
				count++
			}
			if count == 3 {
				return true
			}
		}
	}
	return false
}

func (n *NLS) ToString() string {
	s := ""
	for _, codon := range n.codons {
		s += string(codon.Aa)
	}
	return fmt.Sprintf("%d: %s", n.codons[0].Pos+1, s)
}

func joinNLSes(ns []NLS) string {
	s := make([]string, len(ns))
	for i, n := range ns {
		s[i] = n.ToString()
	}
	return strings.Join(s, ", ")
}

func FindNLS(g *genomes.Genomes, which int) []NLS {
	ret := make([]NLS, 0)
	trans := genomes.Translate(g, which)

	for i := 0; i < len(trans) - 7; i++ {
		nls := NLS{trans[i:i+7]}
		if nls.IsValid() {
			ret = append(ret, nls)
		}
	}
	return ret
}

func LoadShortNames() []string {
	ret := make([]string, 0)
	utils.Lines("../fasta/short_names.txt", func(line string, err error) bool {
		fields := strings.Split(line, " ")
		ret = append(ret, fields[1])
		return true
	})
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	g.Names = LoadShortNames()

	for i := 0; i < g.NumGenomes(); i++ {
		nls := FindNLS(g, i)
		fmt.Printf("%-16s has %d pat7 NLSes at %s\n",
			g.Names[i], len(nls), joinNLSes(nls))
	}
}
