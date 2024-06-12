package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
)

/*
g contains two genomes. Mutate the first one to have total matches with the
second, of which matches need to match. Consider silent mutations only. Returns
the number we found.
*/
func Simulate(g *genomes.Genomes, nd *mutations.NucDistro,
	matches int, total int, its int) int {
	var ret int

iterations:
	for i := 0; i < its; i++ {
		mutant := g.Filter(0)
		mutant.DeepCopy(0)
		mutPositions := mutations.MutateSilent(mutant, nd, total, 1)

		for _, pos := range mutPositions {
			if mutant.Nts[0][pos] != g.Nts[1][pos] {
				continue iterations
			}
		}
		ret++
		fmt.Printf("Found %d/%d %g\n", ret, i+1, float64(ret)/float64(i+1))
	}
	return ret
}
