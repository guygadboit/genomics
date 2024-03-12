package mutations

import (
	"genomics/genomes"
	"genomics/utils"
	"math/rand"
)

// if wantSilent, then make the muts silent, otherwise make them non-silent.
func mutate(genome *genomes.Genomes,
	nucDist *NucDistro, num int, wantSilent bool) int {
	numMuts := 0
	alreadyDone := make(map[int]int)
	nts := genome.Nts[0]

	// Try to mutate silently at pos. Return true if we succeeded.
	tryMutate := func(pos int) bool {
		done, _ := alreadyDone[pos]
		if done != 0 {
			return false
		}

		var env genomes.Environment
		err := env.Init(genome, pos, 1, 0)
		if err != nil {
			// You get an error if pos is not in an ORF
			return false
		}

		existing := nts[pos]
		var replacement byte
		for {
			replacement = nucDist.Random()
			if replacement != existing {
				break
			}
		}

		silent, _ := env.Replace([]byte{replacement})
		success := silent == wantSilent
		if success {
			nts[pos] = replacement
			alreadyDone[pos] = 1
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
func MutateSilent(genome *genomes.Genomes, nucDist *NucDistro, num int) int {
	return mutate(genome, nucDist, num, true)
}

// The same as MutateSilent but make sure the muts are non-silent
func MutateNonSilent(genome *genomes.Genomes,
	nucDist *NucDistro, num int) int {
	return mutate(genome, nucDist, num, false)
}

/*
Returns the number of silent and non-silent mutations in an alignment of
two genomes. Ignores indels.
*/
func CountMutations(g *genomes.Genomes) (int, int) {
	var nonSilent, silent int
	a_nts := g.Nts[0]
	b_nts := g.Nts[1]
	n := g.Length()

	for i := 0; i < n; i++ {
		a := a_nts[i]
		b := b_nts[i]

		if a == b {
			continue
		}

		if a == '-' || b == '-' {
			continue
		}

		err, isSilent, _ := genomes.IsSilent(g, i, 1, 0, 1)
		if err != nil {
			// Ignore anything not in an ORF
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

type Mutation struct {
	Pos    int
	Silent bool
}

func FindMutations(g *genomes.Genomes, a, b int) []Mutation {
	a_nts := g.Nts[a]
	b_nts := g.Nts[b]
	ret := make([]Mutation, 0)

	for i := 0; i < g.Length(); i++ {
		a := a_nts[i]
		b := b_nts[i]

		if a == b {
			continue
		}

		if !utils.IsRegularNt(a) || !utils.IsRegularNt(b) {
			continue
		}

		_, isSilent, _ := genomes.IsSilent(g, i, 1, 0, 1)

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
