package mutations

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"math/rand"
)

/*
genome should contain two genomes. We mutate the first one. If wantSilent, we
make sure the changes are silent relative to the second one. If !wantSilent, we
make sure they are non-silent. Pass 2 into numSeq if you want to mutate two
adjacent nts at once. Returns the positions of the mutations.
*/
func mutate(genome *genomes.Genomes,
	nucDist *NucDistro, num int, numSeq int, wantSilent bool) []int {
	numMuts := 0
	alreadyDone := make(map[int]bool)
	nts := genome.Nts[0]
	maxStart := genome.Length() - numSeq

	// Try to mutate silently (or not silently, as requested) at pos. Return
	// true if we succeeded.
	tryMutate := func(pos int) bool {
		done, _ := alreadyDone[pos]
		if done {
			return false
		}

		existing := nts[pos : pos+numSeq]
		if !utils.IsRegularPattern(existing) {
			return false
		}

		replacement := make([]byte, numSeq)
		other := genome.Nts[1][pos : pos+numSeq]

	findingReplacement:
		for {
			nucDist.RandomSequence(replacement)
			for i := 0; i < len(existing); i++ {
				if other[i] == replacement[i] {
					// It needs to be all different to what's there in the
					// other genome.
					continue findingReplacement
				}
			}
			break
		}

		silent, _, err := genomes.IsSilentWithReplacement(
			genome, pos, 0, 1, replacement)
		if err != nil {
			// You get an error if pos is not in an ORF
			return false
		}

		success := silent == wantSilent
		if success {
			copy(nts[pos:pos+numSeq], replacement)
			for i := pos; i < pos+numSeq; i++ {
				alreadyDone[i] = true
			}
			numMuts++
		}
		return success
	}

mutations:
	for i := 0; i < num; {
		start := rand.Intn(maxStart)

		for j := start; j < maxStart; j++ {
			if tryMutate(j) {
				i++
				continue mutations
			}
		}

		for j := 0; j < start; j++ {
			if tryMutate(j) {
				i++
				continue mutations
			}
		}

		// If we get here it means we've wrapped around the whole genome and
		// not found anywhere to put another silent mutation in. This shouldn't
		// ever happen.
		break
	}
	return utils.FromSet(alreadyDone)
}

// If g has only one genome in it, double it up (so we just mutate it relative
// to itself, which is what the caller intended)
func double(g *genomes.Genomes) *genomes.Genomes {
	if g.NumGenomes() == 1 {
		return g.Filter(0, 0)
	}
	return g
}

/*
Introduce num silent mutations into genome (the first one), selecting nts
randomly from nucDist. Return the number of mutations. If there's only one
genome, just mutate that relative to itself
*/
func MutateSilent(genome *genomes.Genomes,
	nucDist *NucDistro, num int, numSeq int) []int {
	return mutate(double(genome), nucDist, num, numSeq, true)
}

// The same as MutateSilent but make sure the muts are non-silent
func MutateNonSilent(genome *genomes.Genomes,
	nucDist *NucDistro, num int, numSeq int) []int {
	return mutate(double(genome), nucDist, num, numSeq, false)
}

/*
If seqLen is 2, you're looking for two adjacent mutations. Returns the numbers
of silent and non-silent. Ignores indels
*/
func CountSequentialMutations(g *genomes.Genomes, seqLen int) (int, int) {
	var nonSilent, silent int
	a_nts := g.Nts[0]
	b_nts := g.Nts[1]
	n := g.Length()

outer:
	for i := 0; i < n+1-seqLen; i++ {
		a := a_nts[i : i+seqLen]
		b := b_nts[i : i+seqLen]

		for j := 0; j < seqLen; j++ {
			if a[j] == b[j] {
				continue outer
			}
			if a[j] == '-' || b[j] == '-' {
				continue outer
			}
		}

		isSilent, numMuts, _ := genomes.IsSilent(g, i, seqLen, 0, 1)
		if numMuts == 0 {
			continue
		}

		if isSilent {
			silent++
		} else {
			nonSilent++
		}
	}
	return silent, nonSilent
}

/*
Returns the number of silent and non-silent mutations in an alignment of
two genomes. Ignores indels.
*/
func CountMutations(g *genomes.Genomes) (int, int) {
	return CountSequentialMutations(g, 1)
}

type Mutation struct {
	Pos    int
	Silent bool
}

/*
Return a list of the actual mutations between two genomes, rather than just
counting them.
*/
func FindMutations(g *genomes.Genomes, a, b int) []Mutation {
	a_nts := g.Nts[a]
	b_nts := g.Nts[b]
	ret := make([]Mutation, 0)

	for i := 0; i < g.Length(); i++ {
		aNt := a_nts[i]
		bNt := b_nts[i]

		if aNt == bNt {
			continue
		}

		if !utils.IsRegularNt(aNt) || !utils.IsRegularNt(bNt) {
			continue
		}

		isSilent, _, _ := genomes.IsSilent(g, i, 1, a, b)
		ret = append(ret, Mutation{i, isSilent})
	}
	return ret
}

// Wherever a differs from b silently revert a to b. Return the number
// reverted.
func RevertSilent(g *genomes.Genomes, a, b int) int {
	newNts := make([]byte, g.Length())
	var ret int

	for i := 0; i < g.Length(); i++ {
		an := g.Nts[a][i]
		bn := g.Nts[b][i]

		// By default keep what we have. If it's silently different, we will
		// take the b nucleotide below.
		newNts[i] = an

		if an == bn {
			continue
		}

		silent, _, err := genomes.IsSilent(g, i, 1, a, b)
		if err != nil {
			continue
		}

		if silent {
			newNts[i] = bn
			ret++
		}
	}
	g.Nts[a] = newNts
	return ret
}

func Summary(g *genomes.Genomes) {
	n := g.Length()
	for i := 0; i < g.NumGenomes(); i++ {
		g2 := g.Filter(0, i)
		s, ns := CountMutations(g2)
		ss := float64(n-s-ns) / float64(n)
		fmt.Printf("%d. %s| S:%d NS:%d Similarity:%.2f%%\n",
			i, g.Names[i], s, ns, ss*100)
	}
}
