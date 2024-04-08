package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"log"
	"math/rand"
	"os"
	"strings"
)

type Alleles map[byte][]int

func joinInts(ints []int, sep string) string {
	s := make([]string, len(ints))
	for i, v := range ints {
		s[i] = fmt.Sprintf("%d", v)
	}
	return strings.Join(s, sep)
}

// Counts how many "quirks" (differences from the others) in each genome. Maps
// a genome index to a count
type QuirkMap map[int]int

func (q QuirkMap) Combine(other QuirkMap) {
	for k, v := range other {
		q[k] += v
	}
}

// Look for alleles where n viruses share one thing, and everybody else has
// the same other thing.
func (a Alleles) checkNearlyUnique(codon genomes.Codon,
	g *genomes.Genomes, n int) QuirkMap {
	ret := make(QuirkMap)

	if len(a) != 2 {
		return ret
	}

	var us, them byte
	common := make([]int, n)

	for k, v := range a {
		if k == '-' {
			continue
		}
		if len(v) <= n {
			us = k
			copy(common, v)
		} else {
			them = k
		}
	}

	if us != 0 && them != 0 {
		orfs := g.Orfs
		orf, pos, _ := orfs.GetOrfRelative(codon.Pos)

		fmt.Printf("%s:%d: %s got %c, everyone else has %c\n",
			orfs[orf].Name, pos/3+1, joinInts(common, ","), us, them)

		for _, v := range common {
			ret[v] += 1
		}
	}
	return ret
}

// Look for alleles unique to something, but not caring whether the others all
// have the same thing there as each other or are variously different.
func (a Alleles) checkUnique2(codon genomes.Codon,
	g *genomes.Genomes) QuirkMap {
	ret := make(QuirkMap)
	orfs := g.Orfs

	for k, v := range a {
		if k == '-' {
			continue
		}
		if len(v) == 1 {
			orf, pos, _ := orfs.GetOrfRelative(codon.Pos)
			fmt.Printf("%d got %s:%d%c, everyone else something else\n",
				v[0], orfs[orf].Name, pos/3+1, k)
			ret[v[0]] += 1
		}
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
// min SS
func (a Alleles) checkPangolin(codon genomes.Codon,
	g *genomes.Genomes) (int, float64) {
	orfs := g.Orfs

outer:
	for k, v := range a {
		if k == '-' {
			continue
		}
		for _, index := range v {
			if !strings.Contains(g.Names[index], "Pangolin") {
				continue outer
			}
		}
		minSS := minSimilarity(g, v...)
		if len(v) >= 5 {
			orf, pos, _ := orfs.GetOrfRelative(codon.Pos)
			fmt.Printf("%d Pangolins minSS: %.2f got "+
				"%s:%d%c, bats something else\n",
				len(v), minSS, orfs[orf].Name, pos/3+1, k)
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
		if k == '-' {
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
	f, err := os.Create("graph-data.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	for k, v := range qm {
		ss := g.SequenceSimilarity(0, k)
		fmt.Fprintln(w, v, ss)
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

func main() {
	var numSharers int
	var unique bool
	var fasta, orfs string
	var pangolins, controls bool

	flag.IntVar(&numSharers, "n", 1,
		"Number of genomes sharing same unusual thing")
	flag.BoolVar(&unique, "u", false,
		"Just check for unique whatever the others have")
	flag.StringVar(&fasta, "fasta", "../fasta/SARS2-relatives.fasta",
		"Fasta file to use")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs",
		"ORFs file to use")
	flag.BoolVar(&pangolins, "pangolin", false, "Pangolin special")
	flag.BoolVar(&controls, "control", false, "Pangolin controls")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)
	g.RemoveGaps()

	quirks := make(QuirkMap)

	translations := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		translations[i] = genomes.Translate(g, i)
	}

	for i := 0; i < len(translations[0]); i++ {
		alleles := make(Alleles)
		for j := 0; j < g.NumGenomes(); j++ {
			if i >= len(translations[j]) {
				continue
			}
			codon := translations[j][i]
			aa := codon.Aa
			_, there := alleles[aa]
			if !there {
				alleles[aa] = make([]int, 0)
			}
			alleles[aa] = append(alleles[aa], j)
		}

		ref := translations[0][i]
		if pangolins {
			alleles.checkPangolin(ref, g)
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

	if len(quirks) > 0 {
		graphData(quirks, g)
	}
}
