package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

// Show what each genome has in all the places where any of them silently
// differ
func ShowDifferences(g *genomes.Genomes, onlySilent bool) {
outer:
	for i := 0; i < g.Length(); i++ {
		values := make([]byte, g.NumGenomes())
		for j := 0; j < g.NumGenomes(); j++ {
			if onlySilent {
				if j > 0 {
					isSilent, _, _ := genomes.IsSilent(g, i, 1, 0, j)
					if !isSilent {
						continue outer
					}
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
			// Mark places where the civets differ from Tor2
			if values[1] != values[0] {
				fmt.Printf(" interesting")
			}
			fmt.Printf("\n")
		}
	}
}

// Return another genomes object containing the translated spikes
func getSpikes(g *genomes.Genomes) *genomes.Genomes {
	var start, end int
	for _, orf := range g.Orfs {
		if orf.Name == "S" {
			start, end = orf.Start, orf.End
			break
		}
	}

	ret := genomes.NewGenomes(genomes.Orfs{}, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret.Nts[i] = genomes.TranslateAlignedShort(g.Nts[i][start:end])
	}

	return ret
}

func main() {
	g := genomes.LoadGenomes("./ShiSet.fasta", "./SARS1.orfs", false)

	fmt.Println("Nucleotides")
	ShowDifferences(g, true)

	spikes := getSpikes(g)
	fmt.Println("Spikes")
	ShowDifferences(spikes, false)
}
