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
	aNts := g.Nts[a][pos : pos+2]
	bNts := g.Nts[b][pos : pos+2]

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

/*
Return the number of singles, the number of doubles, the number of splittable
doubles, and the counts of which pairs were involved in the doubles
*/
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

/*
If num is 2 return the start of every double mut. If it's 3 return the start of
every triple, etc.
*/
func SequentialMuts(muts []mutations.Mutation,
	num int, requireSilent bool) []mutations.Mutation {
	var muts2 []mutations.Mutation

	// Simplest to get rid of non-silent muts up front
	if requireSilent {
		newMuts := make([]mutations.Mutation, 0)
		for _, mut := range muts {
			if mut.Silent {
				newMuts = append(newMuts, mut)
			}
		}
		muts2 = newMuts
	} else {
		muts2 = muts
	}

	ret := make([]mutations.Mutation, 0)

i:
	for i := 0; i < len(muts2)+1-num; i++ {
		mut := muts2[i]
		for j := 1; j < num; j++ {
			if muts2[i+j].Pos != mut.Pos+j {
				continue i
			}
		}
		ret = append(ret, mut)
	}

	return ret
}

/*
If num is 2, show the doubles, then the singles (that aren't part of doubles)
prefixed by key. If num is 3, start with the triples, etc.
*/
func ShowSequentialMuts(muts []mutations.Mutation, num int,
	requireSilent bool, key string) {
	covered := make(map[int]bool)

	for i := num; i > 0; i-- {
		fmt.Printf("%s %d: ", key, i)
		sMuts := SequentialMuts(muts, i, requireSilent)

		for _, mut := range sMuts {
			if !covered[mut.Pos] {
				fmt.Printf("%d ", mut.Pos)
				// If we found a bunch of doubles first, those positions are
				// "covered", and so don't report them again as singles.
				for j := 0; j < num; j++ {
					covered[mut.Pos+j] = true
				}
			}
		}
		fmt.Printf("\n")
	}
}

func ShowSequentialAll(g *genomes.Genomes, num int, requireSilent bool) {
	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < i; j++ {
			muts := mutations.FindMutations(g, i, j)
			/*
			sim, _ := MakeSimulatedMutant(g, i, j)
			muts = mutations.FindMutations(sim, 0, 1)
			*/
			ShowSequentialMuts(muts, num,
				requireSilent, fmt.Sprintf("%d-%d (%d)", i, j, len(muts)))
		}
	}
}
