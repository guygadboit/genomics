package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"math/rand"
	"sort"
)

type DoubleMap map[string]int

/*
Given a silent double mut starting at pos (so pos and pos+1 are mutated), is it
still silent if you revert either half of it?
*/
func CanSplit(g *genomes.Genomes, a, b int, m mutations.Mutation) bool {
	var env genomes.Environment

	pos := m.Pos
	aNts := g.Nts[a][pos:pos+2]
	bNts := g.Nts[b][pos:pos+2]

	env.Init(g, pos, 2, a)

	// Try replacing just the first nt
	silent, _ := env.Replace([]byte{bNts[0], aNts[1]})
	if silent {
		return true
	}

	// Try replacing just the second
	silent, _ = env.Replace([]byte{aNts[0], bNts[1]})
	if silent {
		return true
	}

	// Replacing one nt at a time doesn't maintain silence. If this happens a
	// lot it implies that the double mutation actually happened as a double
	// mutation
	return false
}

func DoubleMutations(g *genomes.Genomes,
	a, b int, requireSilent bool) (int, int, int, DoubleMap) {
	doubles := make(DoubleMap)
	var doubleCount, singleCount, splittableCount int
	muts := mutations.FindMutations(g, a, b)

	var prev *mutations.Mutation
	for i, mut := range muts {
		if requireSilent && !mut.Silent {
			continue
		}

		if prev != nil && mut.Pos == prev.Pos+1 {
			patterns := []string{
				string(g.Nts[a][prev.Pos : mut.Pos+1]),
				string(g.Nts[b][prev.Pos : mut.Pos+1])}
			sort.Strings(patterns)
			transition := fmt.Sprintf("%s-%s", patterns[0], patterns[1])
			doubles[transition]++
			doubleCount++

			if requireSilent {
				if CanSplit(g, a, b, *prev) {
					splittableCount++
				}
			}


		} else {
			singleCount++
		}
		prev = &muts[i]
	}

	/*
	fmt.Printf("%d-%d %d singles %d doubles of which %d can be split\n",
		a, b, singleCount, doubleCount, splittableCount)
	// fmt.Println(doubles)
	*/
	return singleCount, doubleCount, splittableCount, doubles
}

func (d DoubleMap) Combine(other DoubleMap) {
	for k, v := range other {
		d[k] += v
	}
}

func MontecarloDoubles(g *genomes.Genomes, count int) {
	n := g.NumGenomes()

	dm := make(DoubleMap)
	var totalSingles, totalDoubles, totalSplittable int

	for i := 0; i < count; i++ {
		a, b := rand.Intn(n), rand.Intn(n)
		s, d, sp, m := DoubleMutations(g, a, b, true)
		totalSingles += s
		totalDoubles += d
		totalSplittable += sp
		dm.Combine(m)
	}

	fmt.Printf("%d singles %d doubles (%.2f)\n", totalSingles, totalDoubles,
		float64(totalDoubles)/float64(totalSingles))
	fmt.Printf("%d/%d (%.2f) doubles are splittable\n", totalSplittable,
		totalDoubles, float64(totalSplittable)/float64(totalDoubles))
	fmt.Println(dm)
}

func SimulatePair(g *genomes.Genomes, a, b int) {
	// Get the actual values first
	s, d, _, _ := DoubleMutations(g, a, b, true)
	ratio := float64(d) / float64(s)

	// And now the simulated ones
	g2, silent := MakeSimulatedMutant(g, a, b)

	simS, simD, _, _ := DoubleMutations(g2, 0, 1, true)
	simRatio := float64(simD) / float64(simS)

	fmt.Printf("%d-%d %d %.3f %.3f\n", a, b, silent, ratio, simRatio)
}

// If count == -1 simulate them all exhaustively
func SimulateDoubles(g *genomes.Genomes, count int) {
	n := g.NumGenomes()

	if count == -1 {
		for i := 0; i < n; i++ {
			for j := 0; j < i; j++ {
				SimulatePair(g, i, j)
			}
		}
	} else {
		for i := 0; i < count; i++ {
			a, b := rand.Intn(n), rand.Intn(n)
			SimulatePair(g, a, b)
		}
	}
}
