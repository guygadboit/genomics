package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"math/rand"
	"strings"
)

// The alleles in a given location
type Alleles struct {
	g         *genomes.Genomes
	Pos       int
	Nts       map[string][]int // For the nts in a given codon, which genomes have it?
	Aas       map[byte][]int   // For a given AA, which genomes have it?
	Diversity int              // Number of diff. encodings of majority allele
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

func (a *Alleles) GetOtherNts(index int) map[string]bool {
	ret := make(map[string]bool)
	for k, v := range a.Nts {
		for _, gi := range v {
			if gi != index {
				ret[k] = true
			}
		}
	}
	return ret
}

func (a *Alleles) GetOtherAas(index int) map[byte]bool {
	ret := make(map[byte]bool)
	for k, v := range a.Aas {
		for _, gi := range v {
			if gi != index {
				ret[k] = true
			}
		}
	}
	return ret
}

func (a *Alleles) AASummary(exclude string) string {
	ret := ""
	excl := utils.ToSet([]byte(exclude))
	for k, _ := range a.Aas {
		if !excl[k] {
			ret += string(k)
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

func (a *Alleles) CalcDiversity() {
	var majorityAa byte
	var best int
	for k, v := range a.Aas {
		if len(v) > best {
			best = len(v)
			majorityAa = k
		}
	}

	var diversity int
	for k, _ := range a.Nts {
		if translate(k) == majorityAa {
			diversity++
		}
	}
	a.Diversity = diversity
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
	ret.CalcDiversity()
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
func (a *Alleles) FindUniqueAa(showWhich map[int]bool, cb FoundCB) {
	for k, v := range a.Aas {
		if len(v) == 1 {
			index := v[0]
			others := a.GetOtherAas(index)
			a.HandleMatch(showWhich, index, "UniqueAA", "", k, nil, others)
			cb(index, false)
		}
	}
}

func translate(nts string) byte {
	ret, there := genomes.CodonTable[nts]
	if !there {
		return '-'
	}
	return ret
}

func (a *Alleles) FindSoleOutlierAa(showWhich map[int]bool, cb FoundCB) {
	for k, v := range a.Aas {
		if len(v) == 1 {
			index := v[0]
			others := a.GetOtherAas(index)
			if len(others) == 1 {
				silent := a.HandleMatch(showWhich, index, "SoleOutlierAA",
					"", k, nil, others)
				cb(index, silent)
			}
		}
	}
}

// True if any of the codons in others code for the same AA as codon
func IsSilent(outlier string, others map[string]bool) bool {
	us := translate(outlier)
	for other, _ := range others {
		if us == genomes.CodonTable[other] {
			return true
		}
	}
	return false
}

// Prints it out and returns if it is silent
func (a *Alleles) HandleMatch(showWhich map[int]bool, index int,
	prefix string, outlierNts string, outlierAa byte,
	otherNts map[string]bool, otherAa map[byte]bool) bool {

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
		if us != "" {
			us += fmt.Sprintf(" (%c)", outlierAa)
		} else {
			us = string(outlierAa)
		}
	}

	others := "something else"
	if otherNts != nil && len(otherNts) == 1 {
		others = utils.SetItem(otherNts)
	}
	if otherAa != nil && len(otherAa) == 1 {
		if others != "" {
			others += fmt.Sprintf(" (%c)", utils.SetItem(otherAa))
		} else {
			others = string(utils.SetItem(otherAa))
		}
	}

	var silentString string
	silent := true
	if outlierNts != "" && otherNts != nil && len(otherNts) > 1 {
		if IsSilent(outlierNts, otherNts) {
			silentString = " (silent)"
		} else {
			silentString = " (nonsilent)"
			silent = false
		}
	}

	fmt.Printf("%s: %s:%d: %d got %s, everyone else %s%s diversity: %d\n",
		prefix, orfs[orf].Name, pos/3+1,
		index, us, others, silentString, a.Diversity)

	return silent
}

type FoundCB func(index int, silent bool)

func (a *Alleles) FindUniqueNts(showWhich map[int]bool, cb FoundCB) {
	for nts, v := range a.Nts {
		if len(v) == 1 {
			index := v[0]
			others := a.GetOtherNts(index)
			silent := a.HandleMatch(showWhich,
				index, "UniqueNts", nts, translate(nts), others, nil)
			cb(index, silent)
		}
	}
}

func (a *Alleles) FindSoleOutlierNts(showWhich map[int]bool, cb FoundCB) {
	for nts, v := range a.Nts {
		if len(v) == 1 {
			index := v[0]
			others := a.GetOtherNts(index)
			if len(others) == 1 {
				silent := a.HandleMatch(showWhich, index, "SoleOutlierNts",
					nts, translate(nts), others, nil)
				cb(index, silent)
			}
		}
	}
}

// Update results with the counts based on the alleles at a particular location
func (a *Alleles) Count(g *genomes.Genomes,
	results Results, showWhich map[int]bool) {

	a.FindUniqueAa(showWhich, func(index int, silent bool) {
		results[index].UniqueAas++
	})

	a.FindUniqueNts(showWhich, func(index int, silent bool) {
		if silent {
			results[index].UniqueSilentCodons++
		} else {
			results[index].UniqueNonSilentCodons++
		}
	})

	a.FindSoleOutlierAa(showWhich, func(index int, silent bool) {
		results[index].SoleOutlierAas++
	})

	a.FindSoleOutlierNts(showWhich, func(index int, silent bool) {
		if silent {
			results[index].SoleOutlierSilentCodons++
		} else {
			results[index].SoleOutlierNonSilentCodons++
		}
	})
}

func minSimilarity(g *genomes.Genomes, which ...int) float64 {
	ret := 1.0
	for _, i := range which {
		for _, j := range which {
			if i == j {
				continue
			}
			ss := g.SequenceSimilarity(i, j)
			if ss < ret {
				ret = ss
			}
		}
	}
	return ret
}

/*
Look for alleles that are shared at least minInside times within the group
(which will probably be PCoVs or controls), and not outside it.
*/
func (a *Alleles) CheckGroup(minInside int,
	group map[int]bool, name string) (bool, float64) {
	for k, v := range a.Aas {
		vs := utils.ToSet(v)
		if !utils.IsSubset(vs, group) {
			continue
		}
		if len(v) >= minInside {
			minSS := minSimilarity(a.g, v...)
			orfs := a.g.Orfs
			orf, pos, _ := orfs.GetOrfRelative(a.Pos)
			fmt.Printf("%s:%d %d %s minSS:%.4f got %c, bats "+
				"something else (%s)\n", orfs[orf].Name, pos/3+1,
				len(v), name, minSS, k, a.AASummary(string(k)))
			return true, minSS
		}
	}
	return false, 0
}

func GetPangolins(g *genomes.Genomes, id string) map[int]bool {
	ret := make(map[int]bool)
	for i := 0; i < g.NumGenomes(); i++ {
		if strings.Contains(g.Names[i], id) {
			ret[i] = true
		}
	}
	return ret
}

func GetControls(g *genomes.Genomes, id string) (map[int]bool, int) {
	pangolins := GetPangolins(g, id)
	n := len(pangolins)
	total := g.NumGenomes()
	ret := make(map[int]bool)

	count := 0
	for i := 0; i < n; {
		v := rand.Intn(total)
		if !ret[v] {
			ret[v] = true
			if pangolins[v] {
				count++
			}
			i++
		}
	}

	return ret, count
}

func PangolinControls(alleles *Alleles,
	g *genomes.Genomes, id string, minGroup int) {
	var matches int
	var total, count float64
	for i := 0; i < 1000; i++ {
		controls, pCount := GetControls(g, id)
		if matched, minSS := alleles.CheckGroup(minGroup,
			controls,
			fmt.Sprintf("Controls (%d/%d pangolins)",
				pCount, len(controls))); matched {
			matches++
			total += minSS
			count++
		}
	}
	var meanMinSS float64
	if count != 0.0 {
		meanMinSS = total/count
	}
	fmt.Printf("Matched %d out of 1000 minSS %.2f\n",
		matches, meanMinSS)
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
		pangolinId          string
		minGroup            int
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
	flag.StringVar(&showWhichS, "show", "all", "Only show specified genomes")

	// These are hacks so we can do the Pangolin special on human viruses in
	// the MERS-274 set
	flag.StringVar(&pangolinId, "pango-id", "Pangolin",
		"How to find 'pangolin's")
	flag.IntVar(&minGroup, "min-group", 5, "Min in group")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	var showWhich map[int]bool
	if showWhichS == "all" {
		indices := make([]int, g.NumGenomes())
		for i := 0; i < g.NumGenomes(); i++ {
			indices[i] = i
		}
		showWhich = utils.ToSet(indices)
	} else {
		showWhich = utils.ToSet(utils.ParseInts(showWhichS, ","))
	}

	translations := Translate(g)
	results := make(Results, g.NumGenomes())
	var haveResults bool
	for i, _ := range translations[0] {
		alleles := GetAlleles(g, translations, i)

		if printAlleles {
			alleles.Print()
		}
		if pangolins {
			alleles.CheckGroup(minGroup,
				GetPangolins(g, pangolinId), "Pangolins")
		} else if controls {
			PangolinControls(alleles, g, pangolinId, minGroup)
		} else {
			alleles.Count(g, results, showWhich)
			haveResults = true
		}
	}
	if haveResults {
		SummarizeResults(results)
	}
}
