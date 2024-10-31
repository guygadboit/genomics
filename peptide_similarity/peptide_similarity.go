package main

import (
	"fmt"
	"genomics/genomes"
)

type PeptideSet map[string]bool

func GetPeptides(virus *genomes.Genomes, merLen int) PeptideSet {
	ret := make(PeptideSet)
	for i := 0; i < virus.Length() - merLen; i++ {
		prot := string(virus.Nts[0][i:i+merLen])
		ret[prot] = true
	}
	return ret
}

func Search(peptides PeptideSet, host *genomes.Genomes, merLen int) int {
	ret := 0
	for i := 0; i < host.Length() - merLen; i++ {
		prot := string(host.Nts[0][i:i+merLen])
		if peptides[prot] {
			ret++
		}
	}
	return ret
}

func main() {
	root := "/fs/f/genomes/"

	hosts := []string{
		"human",
		"mouse",
		"chimpanzee",
		"rabbit",
		"pangolin",
	}

	viruses := []string{
		"SARS2",
		"OC43",
		"BANAL-20-52",
		"SARS1",
	}

	for merLen := 6; merLen <= 9; merLen++ {
		for _, v := range viruses {
			vg := genomes.LoadGenomes(
				root+fmt.Sprintf("viruses/%s/%s-prot.fasta.gz", v, v), "", true)
			peptides := GetPeptides(vg, merLen)

			for _, h := range hosts {
				hg := genomes.LoadGenomes(
					root+fmt.Sprintf("%s/%s-prot.fasta.gz", h, h), "", true)
			count := Search(peptides, hg, merLen)
			fmt.Printf("%d %s %s %d\n", merLen, v, h, count)
			}
		}
	}
}
