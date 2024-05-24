package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"math/rand"
)

func OutgroupMontcarlo(its int, numMuts int) {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
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
			for k := 1; k < g.NumGenomes(); k++ {
				if g.Nts[k][pos] == newNt {
					hits++
					break
				}
			}
		}
		total += float64(hits)
		if hits >= 5 {
			passes++
		}
		count++
	}
	fmt.Printf("Mean: %.2f/%d\n", total/count, numMuts)
	fmt.Printf("%d/%d are >= 5 (%.2f)\n",
		passes, its, float64(passes)/float64(its))
}
