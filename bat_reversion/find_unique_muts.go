package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"math/rand"
	"path/filepath"
	"os"
	"strings"
)

// In a given location, for each codon, the list of genomes that have it.
type Alleles map[genomes.Codon][]int

func joinCodons(codons []genomes.Codon, sep string) string {
	s := make([]string, len(codons))
	for i, v := range codons {
		s[i] = v.ToString()
	}
	return strings.Join(s, sep)
}

// True if any of the codons in others code for the same AA as codon
func IsSilent(codon genomes.Codon, others []genomes.Codon) bool {
	for _, other := range others {
		if other.Aa == codon.Aa {
			return true
		}
	}
	return false
}

// Counts how many "quirks" (differences from the others) in each genome. Maps
// a genome index to a count
type QuirkMap struct {
	Silent    map[int]int
	NonSilent map[int]int
}

func (q *QuirkMap) Init() {
	q.Silent = make(map[int]int)
	q.NonSilent = make(map[int]int)
}

func (q QuirkMap) Combine(other QuirkMap) {
	for k, v := range other.Silent {
		q.Silent[k] += v
	}
	for k, v := range other.NonSilent {
		q.NonSilent[k] += v
	}
}

func (q QuirkMap) Summary(names []string) {
	for k, v := range q.Silent {
		name := names[k]
		ns := q.NonSilent[k]
		ratio := float64(ns)/float64(v)
		fmt.Printf("dN/dS: %.2f: %s has %d silent and %d non-silent quirks.\n",
			ratio, name, v, ns)
	}
}

// Look for alleles where n viruses share one thing, and everybody else has
// the same other thing (or a combination of n other things)
func (a Alleles) checkNearlyUnique(codon genomes.Codon,
	g *genomes.Genomes, n int) QuirkMap {
	var ret QuirkMap
	ret.Init()

	var outlier int                             // The index of the outlier
	var us genomes.Codon                        // The allele the outlier has
	alternatives := make([]genomes.Codon, 0, n) // What the others have

	if len(a) > n+1 {
		return ret
	}

	for k, v := range a {
		if k.Aa == '-' {
			continue
		}
		if len(v) == 1 {
			outlier, us = v[0], k
		} else {
			alternatives = append(alternatives, k)
		}
	}

	if us.Aa != 0 && len(alternatives) != 0 {
		orfs := g.Orfs
		orf, pos, _ := orfs.GetOrfRelative(codon.Pos)

		fmt.Printf("%s:%d: %d got %s, everyone else has %s\n",
			orfs[orf].Name, pos/3+1, outlier, us.ToString(),
			joinCodons(alternatives, ","))

		if IsSilent(us, alternatives) {
			ret.Silent[outlier] += 1
		} else {
			ret.NonSilent[outlier] += 1
		}
	}
	return ret
}

// Look for alleles unique to something, but not caring whether the others all
// have the same thing there as each other or are variously different.
func (a Alleles) checkUnique2(codon genomes.Codon,
	g *genomes.Genomes) QuirkMap {
	var ret QuirkMap
	ret.Init()
	orfs := g.Orfs

	for k, v := range a {
		if k.Aa == '-' {
			continue
		}
		if len(v) != 1 {
			continue
		}

		outlier := v[0]

		alternatives := make([]genomes.Codon, 0, 0)
		for other, _ := range a {
			if other != k {
				if other.Aa != '-' {
					alternatives = append(alternatives, other)
				}
			}
		}

		// This might happen if the outlier has an insertion. We don't want to
		// count insertions or deletions here.
		if len(alternatives) == 0 {
			continue
		}

		var silent string

		if IsSilent(k, alternatives) {
			ret.Silent[outlier] += 1
			silent = "(silent)"
		} else {
			ret.NonSilent[outlier] += 1
			silent = "(non-silent)"
		}

		// TODO: just print out what the alternatives are, and how much they
		// differ by.
		orf, pos, _ := orfs.GetOrfRelative(codon.Pos)
		fmt.Printf("%d got %s:%d%s, everyone else something else %s\n",
			outlier, orfs[orf].Name, pos/3+1, k.ToString(), silent)
	}
	return ret
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

// Look for alleles that the pangolin CoVs all share which the bat ones don't.
// Return the number that share the same thing under that criterion, and the
// min SS. FIXME check through this-- haven't done much after porting it to use
// codons.
func (a Alleles) checkPangolin(codon genomes.Codon,
	g *genomes.Genomes, incSC2 bool) (int, float64) {
	orfs := g.Orfs

outer:
	for k, v := range a {
		if k.Aa == '-' {
			continue
		}
		for _, index := range v {
			ok := strings.Contains(g.Names[index], "Pangolin")
			if incSC2 {
				ok = ok || (index == 0)
			}
			if !ok {
				continue outer
			}
		}
		minSS := minSimilarity(g, v...)
		if len(v) >= 5 {
			var sharedWithSC2 string
			if incSC2 {
				if v[0] == 0 {
					sharedWithSC2 = " *and SC2*"
				}
			}
			orf, pos, _ := orfs.GetOrfRelative(codon.Pos)
			fmt.Printf("%d Pangolins minSS: %.2f got "+
				"%s:%d%c, bats something else%s\n",
				len(v), minSS, orfs[orf].Name, pos/3+1, k, sharedWithSC2)

			if sharedWithSC2 != "" {
				for _, index := range v {
					fmt.Printf("%d: %s\n", index, g.Names[index])
				}
			}
		}
		return len(v), minSS
	}
	return 0, 1.0
}

func contains(haystack []int, needle int) bool {
	for _, v := range haystack {
		if v == needle {
			return true
		}
	}
	return false
}

// Look for alleles that 12 CoVs picked at random share which the others don't
func (a Alleles) checkPangolinControl(codon genomes.Codon,
	g *genomes.Genomes, controls []int) (int, float64) {
	// orfs := g.Orfs

outer:
	for k, v := range a {
		if k.Aa == '-' {
			continue
		}
		for _, index := range v {
			if !contains(controls, index) {
				continue outer
			}
		}

		minSS := minSimilarity(g, v...)
		/*
			orf, pos, _ := orfs.GetOrfRelative(codon.Pos)
			fmt.Printf("%d controls minSS: %.2f got %s:%d%c, "+
				"others something else\n",
				len(v), minSS, orfs[orf].Name, pos/3+1, k)
		*/
		return len(v), minSS
	}
	return 0, 1.0
}

func graphData(qm QuirkMap, g *genomes.Genomes) {
	centroid := g.Centroid()

	f, err := os.Create("graph-data.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	// FIXME just doing Silent here
	for k, v := range qm.Silent {
		vec := g.ToVector(k)
		d := utils.VecDistance(vec, centroid)
		fmt.Fprintln(w, v, d, k)
	}

	w.Flush()
	fmt.Println("Wrote graph-data.txt")
}

func randomControls(n int) []int {
	values := make(map[int]bool)
	ret := make([]int, 12)

	for i := 0; i < 12; {
		v := rand.Intn(n)
		if !values[v] {
			ret[i] = v
			i++
			values[v] = true
		}
	}

	return ret
}

func pangolinControls(codon genomes.Codon,
	alleles Alleles, g *genomes.Genomes) {
	var numMatched int
	var total, count float64

	/*
		pangolins := []int{30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41}
		pMatched, pMinSS := alleles.checkPangolinControl(codon, g, pangolins)
		if pMatched > 0 {
			fmt.Printf("Matched %d %f for real pangolins\n", pMatched, pMinSS)
		}
	*/

	for i := 0; i < 1000; i++ {
		matched, minSS := alleles.checkPangolinControl(codon, g,
			randomControls(g.NumGenomes()))
		if matched >= 5 {
			numMatched++
			total += minSS
			count++
		}
	}
	if numMatched > 0 {
		fmt.Printf("Matched %d/1000 minSS %f\n", numMatched, total/count)
	}
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

var SHORT_NAMES []string

func init() {
	SHORT_NAMES = LoadShortNames()
}

func main() {
	var (
		numSharers          int
		unique              bool
		fasta, orfs         string
		pangolins, controls bool
		incSC2              bool
		spikeOnly			bool
	)

	flag.IntVar(&numSharers, "n", 1, "Maximum number of alternatives")
	flag.BoolVar(&unique, "u", false,
		"Just check for unique whatever the others have")
	flag.StringVar(&fasta, "fasta", "../fasta/SARS2-relatives.fasta",
		"Fasta file to use")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs",
		"ORFs file to use")
	flag.BoolVar(&pangolins, "pangolin", false, "Pangolin special")
	flag.BoolVar(&incSC2, "psc2", false, "Include SC2 in Pangolin special")
	flag.BoolVar(&controls, "control", false, "Pangolin controls")
	flag.BoolVar(&spikeOnly, "spike", false, "Spike only")
	flag.Parse()

	useShortNames := filepath.Base(fasta) == "SARS2-relatives.fasta"
	g := genomes.LoadGenomes(fasta, orfs, false)
	g.RemoveGaps()

	var quirks QuirkMap
	quirks.Init()

	translations := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		translations[i] = genomes.Translate(g, i)
	}

	for i := 0; i < len(translations[0]); i++ {
		if spikeOnly {
			t := translations[0]
			if t[i].Pos < 21562 || t[i].Pos >= 25384 {
				continue
			}
		}

		alleles := make(Alleles)
		for j := 0; j < g.NumGenomes(); j++ {
			if i >= len(translations[j]) {
				continue
			}
			codon := translations[j][i]
			/*
				// Uncomment this to treat deletions as not being different alleles
				if codon.Aa == '-' {
					continue
				}
			*/

			_, there := alleles[codon]
			if !there {
				alleles[codon] = make([]int, 0)
			}
			alleles[codon] = append(alleles[codon], j)
		}

		ref := translations[0][i]
		if pangolins {
			alleles.checkPangolin(ref, g, incSC2)
		} else if controls {
			pangolinControls(ref, alleles, g)
		} else if unique {
			q := alleles.checkUnique2(ref, g)
			quirks.Combine(q)
		} else {
			q := alleles.checkNearlyUnique(ref, g, numSharers)
			quirks.Combine(q)
		}
	}

	var names []string
	if useShortNames {
		names = SHORT_NAMES
	} else {
		names = g.Names
	}
	quirks.Summary(names)
	/*
	// FIXME graph both
	if len(quirks.Silent) > 0 {
		graphData(quirks, g)
	}
	*/
}
