package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"math/rand"
	"sort"
)

type DoubleMap map[string]int

func DoubleMutations(g *genomes.Genomes,
	a, b int, requireSilent bool) (int, int, DoubleMap) {
	doubles := make(DoubleMap)
	var doubleCount, singleCount int
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
		} else {
			singleCount++
		}
		prev = &muts[i]
	}

	/*
		fmt.Printf("%d-%d %d singles %d doubles\n", a, b, singleCount, doubleCount)
		fmt.Println(doubles)
	*/
	return singleCount, doubleCount, doubles
}

func (d DoubleMap) Combine(other DoubleMap) {
	for k, v := range other {
		d[k] += v
	}
}

func MontecarloDoubles(g *genomes.Genomes, count int) {
	n := g.NumGenomes()

	dm := make(DoubleMap)
	var totalSingles, totalDoubles int

	for i := 0; i < count; i++ {
		a, b := rand.Intn(n), rand.Intn(n)
		s, d, m := DoubleMutations(g, a, b, true)
		totalSingles += s
		totalDoubles += d
		dm.Combine(m)
	}

	fmt.Printf("%d singles %d doubles (%.2f)\n", totalSingles, totalDoubles,
		float64(totalDoubles)/float64(totalSingles))
	fmt.Println(dm)
}

func SimulatePair(g *genomes.Genomes, a, b int) {
	// Get the actual values first
	s, d, _ := DoubleMutations(g, a, b, true)
	ratio := float64(d) / float64(s)

	// And now the simulated ones
	g2, silent := MakeSimulatedMutant(g, a, b)

	simS, simD, _ := DoubleMutations(g2, 0, 1, true)
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
