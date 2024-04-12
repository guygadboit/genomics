package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

// Show what each genome has in all the places where any of them silently
// differ
func ShowSilent(g *genomes.Genomes) {
outer:
	for i := 0; i < g.Length(); i++ {
		values := make([]byte, g.NumGenomes())
		for j := 0; j < g.NumGenomes(); j++ {
			if j > 0 {
				_, isSilent, _ := genomes.IsSilent(g, i, 1, 0, j)
				if !isSilent {
					continue outer
				}
			}
			values[j] = g.Nts[j][i]
		}
		sValues := utils.ToSet(values)
		if len(sValues) > 1 {
			fmt.Printf("%8d: ", i+1)
			for _, v := range values {
				fmt.Printf("%c", v)
			}
			if values[1] != values[0] {
				fmt.Printf(" interesting")
			}
			fmt.Printf("\n")
		}
	}
}

func main() {
	g := genomes.LoadGenomes("./ShiSet.fasta", "./SARS1.orfs", false)
	ShowSilent(g)
}
