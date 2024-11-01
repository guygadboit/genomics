package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

type PeptideSet map[string]bool

func getPeptides(virus *genomes.Genomes, merLen int) PeptideSet {
	ret := make(PeptideSet)
	for i := 0; i < virus.Length()-merLen; i++ {
		prot := string(virus.Nts[0][i : i+merLen])
		ret[prot] = true
	}
	return ret
}

func Search(peptides PeptideSet, host *genomes.Genomes, merLen int) int {
	matches := make(PeptideSet)
	for i := 0; i < host.Length()-merLen; i++ {
		prot := string(host.Nts[0][i : i+merLen])
		if peptides[prot] {
			matches[prot] = true
		}
	}
	return len(matches)
}

var ROOT string = "/fs/f/genomes/"

func GetPeptides(name string, merLen int, spikeOnly bool) PeptideSet {
	var insert string
	if spikeOnly {
		insert = "S-"
	}
	fname := ROOT + fmt.Sprintf("viruses/%s/%s-%sprot.fasta.gz",
		name, name, insert)

	vg := genomes.LoadGenomes(fname, "", true)
	return getPeptides(vg, merLen)
}

type Library struct {
	hosts   []string
	viruses []string
}

func (l *Library) Init() {
	l.hosts = []string{
		"human",
		"mouse",
		"chimpanzee",
		"rabbit",
		"pangolin",
	}

	l.viruses = []string{
		"SARS2",
		"OC43",
		"BANAL-20-52",
		"RaTG13",
		"SARS1",
	}
}

func (l *Library) Intersections(merLen int, spikeOnly bool) {
	sc2 := GetPeptides(l.viruses[0], merLen, spikeOnly)
	for i, v := range l.viruses {
		if i == 0 {
			continue
		}
		peptides := GetPeptides(v, merLen, spikeOnly)
		ix := len(utils.Intersection(sc2, peptides))
		pct := float64(ix*100) / float64(len(peptides))
		fmt.Printf("SC2 vs %s: %d/%d (%.2f%%) %d-mers in common\n",
			v, ix, len(peptides), pct, merLen)
	}
}

func (l *Library) HostMatches(merLen int, spikeOnly bool) {
	for _, v := range l.viruses {
		peptides := GetPeptides(v, merLen, spikeOnly)

		for _, h := range l.hosts {
			hg := genomes.LoadGenomes(
				ROOT+fmt.Sprintf("%s/%s-prot.fasta.gz", h, h), "", true)
			count := Search(peptides, hg, merLen)

			total := len(peptides)
			pctOutside := (float64(total - count) * 100) / float64(total)
			fmt.Printf("%d %s %s %d/%d %.2f\n", merLen,
				v, h, count, total, pctOutside)
		}
	}
}

func main() {
	var l Library
	l.Init()

	l.Intersections(5, false)
	l.HostMatches(5, false)
}
