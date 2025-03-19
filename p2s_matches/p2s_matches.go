package main

import (
	"fmt"
	"math/rand"
	"genomics/genomes"
)

func Test(g *genomes.Genomes, iterations, length int) int {
	count := 0
	for i := 0; i < iterations; i++ {
		pos := rand.Intn(g.Length()-length)

		matched := true
		for j := 0; j < length; j++ {
			if g.Nts[0][pos+j] != g.Nts[1][pos+j] {
				matched = false
				break
			}
		}
		if matched {
			count++
		}
	}
	return count
}

func main() {
	g := genomes.LoadGenomes("../fasta/WH1-P2S.fasta",
		"../fasta/WH1.orfs", false)

	iterations := 2633 * 100
	matches := Test(g, iterations, 75)
	rate := float64(matches)/float64(iterations)
	fmt.Printf("Matched %d out of %d (%.2f)\n", matches, iterations, rate)
}
