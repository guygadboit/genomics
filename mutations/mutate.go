package mutations

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"math/rand"
)

/*
if wantSilent, then make the muts silent, otherwise make them non-silent
(!wantSilent doesn't mean you don't care, it means you actually want num
non-silent muts) numSeq is 1 if you want single-nt muts. It's 2 if you want
"SDMs" ("sequential double mutations") etc.
*/
func mutate(genome *genomes.Genomes,
	nucDist *NucDistro, num int, numSeq int, wantSilent bool) int {
	numMuts := 0
	alreadyDone := make(map[int]int)
	nts := genome.Nts[0]

	// Try to mutate silently (or not silently, as requested) at pos. Return
	// true if we succeeded.
	tryMutate := func(pos int) bool {
		done, _ := alreadyDone[pos]
		if done != 0 {
			return false
		}

		var env genomes.Environment
		err := env.Init(genome, pos, numSeq, 0)
		if err != nil {
			// You get an error if pos is not in an ORF
			return false
		}

		existing := nts[pos : pos+numSeq]
		replacement := make([]byte, numSeq)

	findingReplacement:
		for {
			nucDist.RandomSequence(replacement)
			for i := 0; i < len(existing); i++ {
				if existing[i] == replacement[i] {
					// It needs to be all different
					continue findingReplacement
				}
			}
			break
		}

		silent, _ := env.Replace(replacement)
		success := silent == wantSilent
		if success {
			copy(nts[pos:pos+numSeq], replacement)
			for i := pos; i < pos+numSeq; i++ {
				alreadyDone[i] = 1
			}
			numMuts++
		}
		return success
	}

mutations:
	for i := 0; i < num; {
		start := rand.Intn(genome.Length())

		for j := start; j < genome.Length(); j++ {
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
	return numMuts
}

/*
Introduce num silent mutations into genome (the first one), selecting nts
randomly from nucDist. Return the number of mutations
*/
func MutateSilent(genome *genomes.Genomes,
	nucDist *NucDistro, num int, numSeq int) int {
	return mutate(genome, nucDist, num, numSeq, true)
}

// The same as MutateSilent but make sure the muts are non-silent
func MutateNonSilent(genome *genomes.Genomes,
	nucDist *NucDistro, num int, numSeq int) int {
	return mutate(genome, nucDist, num, numSeq, false)
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

		_, isSilent, _ := genomes.IsSilent(g, i, seqLen, 0, 1)

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

		_, isSilent, _ := genomes.IsSilent(g, i, 1, a, b)
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

		err, silent, _ := genomes.IsSilent(g, i, 1, a, b)
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
		fmt.Printf("<%s> S:%d NS:%d Similarity: %.2f%%\n",
			g.Names[i], s, ns, ss*100)
	}
}
