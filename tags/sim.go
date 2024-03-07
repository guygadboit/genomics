package main

import (
	"genomics/genomes"
	"genomics/mutations"
)

/*
Given the two genomes a, b in g, count how many silent muts from a to b, and
return a new alignment of c and a. c is a new simulated mutant with the same
number of silent muts, and a is the same a (just shallow copied-- we aren't
modifying it).
*/
func MakeSimulatedMutant(g *genomes.Genomes, a, b int) *genomes.Genomes {
	ret := genomes.NewGenomes(g.Orfs, 2)
	n := g.Length()

	// Copy the first genome in-- we are going to mutate it
	ret.Nts[0] = make([]byte, n)
	copy(ret.Nts[0], g.Nts[a])

	// The second one can be a shallow copy of b while we count the number of
	// muts and get the nt distribution of both of them.
	ret.Nts[1] = g.Nts[b]
	silent, _ := mutations.CountMutations(g)
	nd := mutations.NewNucDistro(ret)

	mutations.MutateSilent(ret, nd, silent)

	// Now put the original back in the second slot because we're going to look
	// for tags between the two.
	ret.Nts[1] = g.Nts[a]

	ret.Names[0] = "Simulated Mutant"
	ret.Names[1] = g.Names[0]

	return ret
}
