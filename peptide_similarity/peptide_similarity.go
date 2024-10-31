package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

type PeptideSet map[string]bool

func GetPeptides(virus *genomes.Genomes, merLen int) PeptideSet {
	ret := make(PeptideSet)
	for i := 0; i < virus.Length()-merLen; i++ {
		prot := string(virus.Nts[0][i : i+merLen])
		ret[prot] = true
	}
	return ret
}

func Search(peptides PeptideSet, host *genomes.Genomes, merLen int) int {
	ret := 0
	for i := 0; i < host.Length()-merLen; i++ {
		prot := string(host.Nts[0][i : i+merLen])
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
		"RaTG13",
		"SARS1",
	}

	getPeptides := func(name string, merLen int) PeptideSet {
		vg := genomes.LoadGenomes(
			root+fmt.Sprintf("viruses/%s/%s-prot.fasta.gz",
				name, name), "", true)
		return GetPeptides(vg, merLen)
	}

	for merLen := 5; merLen <= 9; merLen++ {
		sc2 := getPeptides("SARS2", merLen)
		for i, v := range viruses {
			if i == 0 {
				continue
			}
			peptides := getPeptides(v, merLen)
			ix := len(utils.Intersection(sc2, peptides))
			pct := float64(ix*100) / float64(len(peptides))
			fmt.Printf("SC2 vs %s: %d/%d (%.2f%%) %d-mers in common\n",
				v, ix, len(peptides), pct, merLen)
		}
	}
	return

	for merLen := 6; merLen <= 9; merLen++ {
		for _, v := range viruses {
			peptides := getPeptides(v, merLen)

			for _, h := range hosts {
				hg := genomes.LoadGenomes(
					root+fmt.Sprintf("%s/%s-prot.fasta.gz", h, h), "", true)
				count := Search(peptides, hg, merLen)
				fmt.Printf("%d %s %s %d\n", merLen, v, h, count)
			}
		}
	}
}
