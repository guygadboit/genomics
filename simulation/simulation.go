package simulation

import (
	"genomics/genomes"
	"genomics/mutations"
)

/*
Given the two genomes a, b in g, count how many silent muts from a to b, and
return a new alignment of c and a. c is a new simulated mutant with the same
number of silent muts, and a is the same a (just shallow copied-- we aren't
modifying it). Also returns the number of silent muts
*/
func MakeSimulatedMutant(g *genomes.Genomes, a, b int) (*genomes.Genomes, int) {
	// First obtain a nucleotide distribution based on both of them
	g2 := g.Filter(a, b)
	silent, _ := mutations.CountMutations(g2)
	nd := mutations.NewNucDistro(g2)

	ret := g.Filter(a, b)

	// Now take out the silent muts
	ret.DeepCopy(0)	// Take a copy so we don't mess up g
	mutations.RevertSilent(ret, 0, 1)

	// And put them back randomly. We're interested to see if this gives
	// different results.
	mutations.MutateSilent(ret, nd, 1, silent)

	ret.Names[0] = "Simulated Mutant"
	return ret, silent
}

/*
Count the number of silent and non-silent muts between a, b, and return
something that contains a mutated version of a, with the same numbers of each,
but distributed evenly, and the original a.
*/
func MakeSimulatedMutant2(g *genomes.Genomes,
	a, b int) (*genomes.Genomes, int) {
	g2 := g.Filter(a, b)
	silent, nonSilent := mutations.CountMutations(g2)
	nd := mutations.NewNucDistro(g2)

	ret := g.Filter(a, a)
	ret.DeepCopy(0)
	mutations.MutateSilent(ret, nd, 1, silent)
	mutations.MutateNonSilent(ret, nd, 1, nonSilent)

	ret.Names[0] = "Simulated Mutant"
	return ret, silent
}
