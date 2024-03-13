package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/simulation"
)

/*
Represents a concentration of NumMuts mutations in a sequence of nts Length
long. So if Length and NumMuts are 2,2 this looks for "doubles" (mutations next
to each other). If they're 6,4 it looks for what you were calling "tags"
*/
type Concentration struct {
	Pos     int
	Length  int
	NumMuts int
	Silent  bool
}

func FindConcentrations(g *genomes.Genomes, length int, minMuts int,
	requireSilent bool) []Concentration {
	ret := make([]Concentration, 0)

positions:
	for i := 0; i < g.Length()-length; i++ {
		for j := 0; j < g.NumGenomes(); j++ {
			for k := 0; k < j; k++ {
				err, silent, numMuts := genomes.IsSilent(g, i, length, j, k)
				if err != nil {
					continue
				}
				ok := silent || !requireSilent
				if ok && numMuts >= minMuts {
					// fmt.Printf("Found pattern at %d\n", i)
					ret = append(ret,
						Concentration{i, length, numMuts, requireSilent})
					continue positions
				}
			}
		}
	}
	return ret
}

func CreateHighlights(concentrations []Concentration) []genomes.Highlight {
	ret := make([]genomes.Highlight, len(concentrations))
	for i, p := range concentrations {
		ret[i] = genomes.Highlight{p.Pos, p.Pos + p.Length, 'v'}
	}
	return ret
}

// Compare counts of concentrations to simulations
func CompareToSim(g *genomes.Genomes, length int, minMuts int,
	requireSilent bool) {
	n := g.NumGenomes()
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			g2 := g.Filter(i, j)
			concs := FindConcentrations(g2, length, minMuts, requireSilent)
			realCount := len(concs)

			simG, numSilent := simulation.MakeSimulatedMutant(g2, 0, 1)
			concs = FindConcentrations(simG, length, minMuts, requireSilent)
			simCount := len(concs)

			fmt.Printf("%d.%d %d-%d: %d %d %.2f (%d)\n",
			length, minMuts, i, j, realCount, simCount,
			float64(simCount)/float64(realCount), numSilent)
		}
	}
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)

	CompareToSim(g, 3, 3, true)
}
