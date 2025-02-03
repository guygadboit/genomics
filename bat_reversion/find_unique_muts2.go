package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
)

// The alleles in a given location
type Alleles struct {
	g   *genomes.Genomes
	Pos int
	Nts map[string][]int // For the nts in a given codon, which genomes have it?
	Aas map[byte][]int   // For a given AA, which genomes have it?
}

func (a *Alleles) Print() {
	orfs := a.g.Orfs
	orf, pos, _ := orfs.GetOrfRelative(a.Pos)

	fmt.Printf("%s:%d: ", orfs[orf].Name, pos/3+1)
	for k, v := range a.Nts {
		fmt.Printf("%s: ", k)
		fmt.Println(v)
	}
	for k, v := range a.Aas {
		fmt.Printf("%c ", k)
		fmt.Println(v)
	}
	fmt.Printf("\n")
}

func NewAlleles(g *genomes.Genomes, pos int) *Alleles {
	var ret Alleles
	ret.g = g
	ret.Pos = pos
	ret.Nts = make(map[string][]int)
	ret.Aas = make(map[byte][]int)
	return &ret
}

func (a *Alleles) Add(codon *genomes.Codon, genomeIndex int) {
	if codon.Aa == '-' {
		return
	}
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

func GetAlleles(g *genomes.Genomes,
	translations []genomes.Translation, codonIndex int) *Alleles {
	ret := NewAlleles(g, translations[0][codonIndex].Pos)
	for i, _ := range translations {
		codon := translations[i][codonIndex]
		ret.Add(&codon, i)
	}
	return ret
}

type Result struct {
	UniqueAas                  int
	SoleOutlierAas             int
	UniqueSilentCodons         int
	UniqueNonSilentCodons      int
	SoleOutlierSilentCodons    int
	SoleOutlierNonSilentCodons int
}

// One fo each genome
type Results []Result

// If there is a unique Aa in there, return the genome that has it. Otherwise
// -1
func (a *Alleles) FindUniqueAa() int {
	index := -1
	var us byte

	for k, v := range a.Aas {
		if len(v) == 1 {
			index = v[0]
			us = k
			break
		}
	}
	if index == -1 {
		return -1
	}

	orfs := a.g.Orfs
	orf, pos, _ := orfs.GetOrfRelative(a.Pos)
	fmt.Printf("Unique Aa: %s:%d: %d got %c, "+
		"everyone else something else.\n", orfs[orf].Name, pos/3+1, index, us)
	return index
}

func translate(nts string) byte {
	ret, there := genomes.CodonTable[nts]
	if !there {
		return '-'
	}
	return ret
}

func silentString(silent bool) string {
	if silent {
		return "silent"
	} else {
		return "nonsilent"
	}
}

func (a *Alleles) FindSoleOutlierAa() int {
	others := make(map[byte]bool)
	var us, alt byte
	index := -1

	for k, v := range a.Aas {
		if len(v) == 1 {
			index = v[0]
			us = k
		} else {
			others[k] = true
			alt = k
			if len(others) > 1 {
				return -1
			}
		}
	}

	if index == -1 {
		return index
	}

	orfs := a.g.Orfs
	orf, pos, _ := orfs.GetOrfRelative(a.Pos)
	fmt.Printf("SoleOutlier Aa: %s:%d: %d got %c, "+
		"everyone else %c.\n", orfs[orf].Name, pos/3+1, index, us, alt)

	return index
}

// True if any of the codons in others code for the same AA as codon
func IsSilent(outlier string, others []string) (bool, byte) {
	us := translate(outlier)
	for _, other := range others {
		if us == genomes.CodonTable[other] {
			return true, us
		}
	}
	return false, us
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

	silent, aa := IsSilent(outlier, others)

	orfs := a.g.Orfs
	orf, pos, _ := orfs.GetOrfRelative(a.Pos)
	fmt.Printf("Unique Nts: %s:%d: %d got %s (%c) (%s), "+
		"everyone else something else.\n",
		orfs[orf].Name, pos/3+1, index, outlier, aa, silentString(silent))

	return index, silent
}

func (a *Alleles) FindSoleOutlierNts() (int, bool) {
	others := make(map[string]bool)
	index := -1
	var outlier, alt string

	for nts, v := range a.Nts {
		if len(v) == 1 {
			index = v[0]
			outlier = nts
		} else {
			others[nts] = true
			if len(others) > 1 {
				return -1, false
			}
			alt = nts
		}
	}

	if index == -1 {
		return -1, false
	}

	aa := translate(outlier)
	altAa := translate(alt)
	silent := aa == altAa

	orfs := a.g.Orfs
	orf, pos, _ := orfs.GetOrfRelative(a.Pos)
	fmt.Printf("SoleOutlier Nts: %s:%d: %d got %s (%c), "+
		"everyone else %s (%c) (%s).\n",
		orfs[orf].Name, pos/3+1, index, outlier, aa, alt, altAa, silentString(silent))

	return index, silent
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

	index, silent = a.FindSoleOutlierNts()
	if index != -1 {
		if silent {
			results[index].SoleOutlierSilentCodons++
		} else {
			results[index].SoleOutlierNonSilentCodons++
		}
	}
}

func main() {
	var (
		unique              bool
		fasta, orfs         string
		pangolins, controls bool
		spikeOnly           bool
		exclude             string
		printAlleles		bool
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
	flag.BoolVar(&printAlleles, "print", false, "Print all alleles")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)
	translations := Translate(g)
	results := make(Results, g.NumGenomes())
	for i, _ := range translations[0] {
		alleles := GetAlleles(g, translations, i)
		if printAlleles {
			alleles.Print()
		}
		alleles.Count(g, results)
	}

	// TODO: Print a summary of counts at the end.
}
