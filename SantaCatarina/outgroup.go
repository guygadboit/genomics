package main

import (
	"fmt"
	"genomics/mutations"
	"genomics/genomes"
	"genomics/database"
	"math/rand"
)

func OutgroupMontcarlo(g *genomes.Genomes, its int, numMuts int) {
	nd := mutations.NewNucDistro(g)

	total, count, passes := 0.0, 0.0, 0

	for i := 0; i < its; i++ {
		hits := 0
		for j := 0; j < numMuts; j++ {
			pos := rand.Intn(g.Length())
			existing := g.Nts[0][pos]

			var newNt byte
			for {
				newNt = nd.Random()
				if newNt != existing {
					break
				}
			}
			// Just looking at some of the PCoVs here
			good := true
			for k := 35; k < 42; k++ {
				if g.Nts[k][pos] != newNt {
					good = false
					break
				}
			}
			if good {
				hits++
			}
		}
		total += float64(hits)
		if hits >= 3 {
			passes++
		}
		count++
	}
	fmt.Printf("Mean: %.2f/%d\n", total/count, numMuts)
	fmt.Printf("%d/%d are >= 3 (%.2f)\n",
		passes, its, float64(passes)/float64(its))
}

func ShowOutgroupMatches(g *genomes.Genomes) {
	muts := database.ParseMutations("T5929G,T8601C,A8651C,G16206A,T19218G")
	// muts := database.ParseMutations("T5929G,A8651C,G16206A")

	for i := 42; i < 44; i++ {
	// for i := 1; i < g.NumGenomes(); i++ {
		count := 0
		for _, m := range muts {
			if g.Nts[i][m.Pos-1] == m.To {
				fmt.Printf("Have %s\n", m.ToString())
				count++
			}
		}
		if count > 0 {
			fmt.Printf("%d %s has %d\n", i, g.Names[i], count)
		}
	}
}
