package simulation

import (
	"genomics/genomes"
	"genomics/mutations"
)

type MutantFunc func(*genomes.Genomes,
	int, int, *mutations.NucDistro) (*genomes.Genomes, int)

/*
Given the two genomes a, b in g, count how many silent muts from a to b, and
return a new alignment of c and a. c is a new simulated mutant with the same
number of silent muts, and a is the same a (just shallow copied-- we aren't
modifying it). Also returns the number of silent muts. nd can be nil which
means we make our own.
*/
func MakeSimulatedMutant(g *genomes.Genomes,
	a, b int, nd *mutations.NucDistro) (*genomes.Genomes, int) {
	// First obtain a nucleotide distribution based on both of them
	g2 := g.Filter(a, b)
	silent, _ := mutations.CountMutations(g2)

	if nd == nil {
		nd = mutations.NewNucDistro(g2)
	}

	ret := g.Filter(a, b)

	// Now take out the silent muts
	ret.DeepCopy(0) // Take a copy so we don't mess up g
	mutations.RevertSilent(ret, 0, 1)

	// And put them back randomly. We're interested to see if this gives
	// different results.
	mutations.MutateSilent(ret, nd, silent, 1)

	ret.Names[0] = "Simulated Mutant"
	return ret, silent
}

/*
Count the number of silent and non-silent muts between a, b, and return
something that contains a mutated version of a, with the same numbers of each,
but distributed evenly, and the original a.
*/
func MakeSimulatedMutant2(g *genomes.Genomes,
	a, b int, nd *mutations.NucDistro) (*genomes.Genomes, int) {
	g2 := g.Filter(a, b)
	silent, nonSilent := mutations.CountMutations(g2)

	if nd == nil {
		nd = mutations.NewNucDistro(g2)
	}

	ret := g.Filter(a, a)
	ret.DeepCopy(0)
	mutations.MutateSilent(ret, nd, silent, 1)
	mutations.MutateNonSilent(ret, nd, nonSilent, 1)

	ret.Names[0] = "Type 2 Simulated Mutant"
	return ret, silent
}

/*
Redistribute the silent mutations, but find the doubles first, and then put
them back (evenly) and then the remaining singles.
*/
func MakeSimulatedMutant3(g *genomes.Genomes,
	a, b int, nd *mutations.NucDistro) (*genomes.Genomes, int) {
	g2 := g.Filter(a, b)
	sDoubles, _ := mutations.CountSequentialMutations(g2, 2)
	sSingles, _ := mutations.CountMutations(g2)

	if nd == nil {
		nd = mutations.NewNucDistro(g2)
	}

	ret := g.Filter(a, b)

	// Now take out all the silent muts
	ret.DeepCopy(0)
	mutations.RevertSilent(ret, 0, 1)

	// Now put back in the right number of doubles
	mutations.MutateSilent(ret, nd, sDoubles, 2)

	// And then any remaining singles
	mutations.MutateSilent(ret, nd, sSingles-sDoubles*2, 1)

	ret.Names[0] = "Type 3 Simulated Mutant"
	return ret, sSingles
}

func MakeSimulatedMutant4(g *genomes.Genomes,
	a, b int, nd *mutations.NucDistro) (*genomes.Genomes, int) {
	g2 := g.Filter(a, b)
	sDoubles, nsDoubles := mutations.CountSequentialMutations(g2, 2)
	sSingles, nsSingles := mutations.CountMutations(g2)

	if nd == nil {
		nd = mutations.NewNucDistro(g2)
	}

	ret := g.Filter(a, a)
	ret.DeepCopy(0)

	mutations.MutateSilent(ret, nd, sDoubles, 2)
	mutations.MutateNonSilent(ret, nd, nsDoubles, 1)

	mutations.MutateSilent(ret, nd, sSingles, 1)
	mutations.MutateNonSilent(ret, nd, nsSingles, 1)

	ret.Names[0] = "Type 4 Simulated Mutant"
	return ret, sSingles
}
