package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
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

func (a *Alleles) GetOtherNts(index int) []string {
	ret := make([]string, 0)
	for k, v := range a.Nts {
		for _, gi := range v {
			if gi != index {
				ret = append(ret, k)
			}
		}
	}
	return ret
}

func (a *Alleles) GetOtherAas(index int) []byte {
	ret := make([]byte, 0)
	for k, v := range a.Aas {
		for _, gi := range v {
			if gi != index {
				ret = append(ret, k)
			}
		}
	}
	return ret
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

func SummarizeResults(results Results) {
	fmt.Printf("Index UniqueAas SoleOutlierAas UniqueSilentCodons " +
		"UniqueNonSilentCodons SoleOutlierSilentCodons " +
		"SoleOutlierNonSilentCodons dN/dS Codons\n")
	for i, result := range results {
		ratio := float64(
			result.UniqueNonSilentCodons) / float64(result.UniqueSilentCodons)

		fmt.Printf("%d: %d %d %d %d %d %d %.3f\n", i,
			result.UniqueAas,
			result.SoleOutlierAas,
			result.UniqueSilentCodons,
			result.UniqueNonSilentCodons,
			result.SoleOutlierSilentCodons,
			result.SoleOutlierNonSilentCodons,
			ratio)
	}
}

// If there is a unique Aa in there, return the genome that has it. Otherwise
// -1
func (a *Alleles) FindUniqueAa(showWhich map[int]bool) int {
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

	if showWhich[index] {
		orfs := a.g.Orfs
		orf, pos, _ := orfs.GetOrfRelative(a.Pos)
		fmt.Printf("Unique Aa: %s:%d: %d got %c, "+
			"everyone else something else.\n", orfs[orf].Name, pos/3+1, index, us)
	}
	return index
}

func translate(nts string) byte {
	ret, there := genomes.CodonTable[nts]
	if !there {
		return '-'
	}
	return ret
}

func (a *Alleles) FindSoleOutlierAa(showWhich map[int]bool) int {
	others := make(map[byte]bool)
	var us byte
	index := -1
	var alt byte = '-'

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

	if showWhich[index] {
		orfs := a.g.Orfs
		orf, pos, _ := orfs.GetOrfRelative(a.Pos)
		fmt.Printf("SoleOutlier Aa: %s:%d: %d got %c, "+
			"everyone else %c.\n", orfs[orf].Name, pos/3+1, index, us, alt)
	}

	return index
}

// True if any of the codons in others code for the same AA as codon
func IsSilent(outlier string, others []string) bool {
	us := translate(outlier)
	for _, other := range others {
		if us == genomes.CodonTable[other] {
			return true
		}
	}
	return false
}

// Prints it out and returns if it is silent
func (a *Alleles) HandleMatch(showWhich map[int]bool, index int,
	prefix string, outlierNts string, outlierAa byte,
	otherNts []string, otherAa []byte) bool {

	if !showWhich[index] {
		return true
	}

	orfs := a.g.Orfs
	orf, pos, _ := orfs.GetOrfRelative(a.Pos)

	var us string
	if outlierNts != "" {
		us = outlierNts
	}
	if outlierAa != 0 {
		us += fmt.Sprintf(" (%c)", outlierAa)
	}

	others := "something else"
	if otherNts != nil && len(otherNts) == 1 {
		others = otherNts[0]
	}
	if otherAa != nil && len(otherAa) > 0 {
		others += fmt.Sprintf(" (%c)", otherAa[0])
	}

	var silentString string
	silent := true
	if outlierNts != "" && otherNts != nil && len(otherNts) > 1 {
		if IsSilent(outlierNts, otherNts) {
			silentString = "silent"
		} else {
			silentString = "nonsilent"
			silent = false
		}
	}

	fmt.Printf("%s: %s:%d: %d got %s, everyone else %s (%s).\n", prefix,
		orfs[orf].Name, pos/3+1, index, us, others, silentString)

	return silent
}


type FoundCB func(index int, silent bool)

func (a *Alleles) FindUniqueNts(showWhich map[int]bool, cb FoundCB) {
	index := -1

	for nts, v := range a.Nts {
		if len(v) == 1 {
			index = v[0]
			others := a.GetOtherNts(index)
			silent := a.HandleMatch(showWhich, 
				index, "UniqueNts", nts, translate(nts), others, nil)
			cb(index, silent)
		}
	}
}

/*
func (a *Alleles) FindSoleOutlierNts(showWhich map[int]bool) (int, bool) {
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

	if showWhich[index] {
		orfs := a.g.Orfs
		orf, pos, _ := orfs.GetOrfRelative(a.Pos)
		fmt.Printf("SoleOutlier Nts: %s:%d: %d got %s (%c), "+
			"everyone else %s (%c) (%s).\n",
			orfs[orf].Name, pos/3+1, index, outlier,
			aa, alt, altAa, silentString(silent))
	}

	return index, silent
}
*/

// Update results with the counts based on the alleles at a particular location
func (a *Alleles) Count(g *genomes.Genomes,
	results Results, showWhich map[int]bool) {

		/*
	index := a.FindUniqueAa(showWhich)
	if index != -1 {
		results[index].UniqueAas++
	}

	index = a.FindSoleOutlierAa(showWhich)
	if index != -1 {
		results[index].SoleOutlierAas++
	}
	*/

	a.FindUniqueNts(showWhich, func(index int, silent bool) {
		if silent {
			results[index].UniqueSilentCodons++
		} else {
			results[index].UniqueNonSilentCodons++
		}
	})

	/*

	index, silent = a.FindSoleOutlierNts(showWhich)
	if index != -1 {
		if silent {
			results[index].SoleOutlierSilentCodons++
		} else {
			results[index].SoleOutlierNonSilentCodons++
		}
	}
	*/
}

func main() {
	var (
		unique              bool
		fasta, orfs         string
		pangolins, controls bool
		spikeOnly           bool
		exclude             string
		printAlleles        bool
		showWhichS          string
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
	flag.StringVar(&showWhichS, "show", "0", "Only show specified genomes")
	flag.Parse()

	showWhich := utils.ToSet(utils.ParseInts(showWhichS, ","))

	g := genomes.LoadGenomes(fasta, orfs, false)
	translations := Translate(g)
	results := make(Results, g.NumGenomes())
	for i, _ := range translations[0] {
		alleles := GetAlleles(g, translations, i)
		if printAlleles {
			alleles.Print()
		}
		alleles.Count(g, results, showWhich)
	}
	SummarizeResults(results)
}
