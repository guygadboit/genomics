package main

import (
	"flag"
	"genomics/genomes"
)

// The alleles in a given location
type Alleles struct {
	Pos int
	Nts map[string][]int // For the nts in a given codon, which genomes have it?
	Aas map[byte][]int   // For a given AA, which genomes have it?
}

func NewAlleles(pos int) *Alleles {
	var ret Alleles
	ret.Pos = pos
	ret.Nts = make(map[string][]int)
	ret.Aas = make(map[byte][]int)
	return &ret
}

func (a *Alleles) Add(codon *genomes.Codon, genomeIndex int) {
	_, there := a.Nts[codon.Nts]
	if !there {
		a.Nts[codon.Nts] = make([]int, 0)
	}
	a.Nts[codon.Nts] = append(a.Nts[codon.Nts], genomeIndex)

	_, there = a.Aas[codon.Aa]
	if !there {
		a.Aas[codon.Aa] = make([]int, 0)
	}
	a.Aas[codon.Aa] = append(a.Aas[codon.Aa], genomeIndex)
}

func Translate(g *genomes.Genomes) []genomes.Translation {
	ret := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret[i] = genomes.Translate(g, i)
	}
	return ret
}

func GetAlleles(translations []genomes.Translation, codonPos int) *Alleles {
	ret := NewAlleles(codonPos)
	for i, _ := range translations {
		codon := translations[i][codonPos]
		ret.Add(&codon, i)
	}
	return ret
}

type Result struct {
	UniqueAas             int
	SoleOutlierAas        int
	UniqueSilentCodons    int
	UniqueNonSilentCodons int
	SoleOutlierCodons     int
}

// One fo each genome
type Results []Result

// If there is a unique Aa in there, return the genome that has it. Otherwise
// -1
func (a *Alleles) FindUniqueAa() int {
	for _, v := range a.Aas {
		if len(v) == 1 {
			return v[0]
		}
	}
	return -1
}

func (a *Alleles) FindSoleOutlierAa() int {
	others := make(map[byte]bool)
	ret := -1

	for k, v := range a.Aas {
		if len(v) == 1 {
			ret = v[0]
		} else {
			others[k] = true
			if len(others) > 1 {
				return -1
			}
		}
	}

	return ret
}

// True if any of the codons in others code for the same AA as codon
func IsSilent(outlier string, others []string) bool {
	for _, other := range others {
		if genomes.CodonTable[outlier] == genomes.CodonTable[other] {
			return true
		}
	}
	return false
}

func (a *Alleles) FindUniqueNts() (int, bool) {
	index := -1

	var outlier string
	others := make([]string, 0)

	for nts, v := range a.Nts {
		if len(v) == 1 {
			index = v[0]
			outlier = nts
			break
		} else {
			others = append(others, nts)
		}
	}
	if index == -1 {
		return -1, false
	}
	return index, IsSilent(outlier, others)
}

func (a *Alleles) FindSoleOutlierNts() int {
	others := make(map[string]bool)
	ret := -1

	for k, v := range a.Nts {
		if len(v) == 1 {
			ret = v[0]
		} else {
			others[k] = true
			if len(others) > 1 {
				return -1
			}
		}
	}

	return ret
}

// Update results with the counts based on the alleles at a particular location
func (a *Alleles) Count(g *genomes.Genomes, results Results) {

	index := a.FindUniqueAa()
	if index != -1 {
		results[index].UniqueAas++
	}

	index = a.FindSoleOutlierAa()
	if index != -1 {
		results[index].SoleOutlierAas++
	}

	index, silent := a.FindUniqueNts()
	if index != -1 {
		if silent {
			results[index].UniqueSilentCodons++
		} else {
			results[index].UniqueNonSilentCodons++
		}
	}

	index = a.FindSoleOutlierNts()
	if index != -1 {
		results[index].SoleOutlierCodons++
	}
}

func main() {
	var (
		unique              bool
		fasta, orfs         string
		pangolins, controls bool
		spikeOnly           bool
		exclude             string
	)

	flag.BoolVar(&unique, "u", false,
		"Just check for unique whatever the others have")
	flag.StringVar(&fasta, "fasta", "../fasta/SARS2-relatives.fasta",
		"Fasta file to use")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs",
		"ORFs file to use")
	flag.BoolVar(&pangolins, "pangolin", false, "Pangolin special")
	flag.BoolVar(&controls, "control", false, "Pangolin controls")
	flag.BoolVar(&spikeOnly, "spike", false, "Spike only")
	flag.StringVar(&exclude, "exclude", "", "Indices to exclude")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)
	translations := Translate(g)
	var results Results
	for i, _ := range translations {
		alleles := GetAlleles(translations, i)
		alleles.Count(g, results)
	}
}
