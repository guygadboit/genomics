package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"math"
)

const SITE_LEN = 6

// Represents anywhere you might end up with a site because somebody in the set
// has one.
type Locus struct {
	Pos   int
	P     float64 // Probability of SC2 having a site here after mutations
}

func (l *Locus) Init(g *genomes.Genomes, pos int, window int) {
	l.Pos = pos

	// Now the harder part: compute P
	for i := 1; i < g.NumGenomes(); i++ {
		start := max(0, pos-window)
		end := min(g.Length(), pos+window)
		possibleMuts := mutations.PossibleSilentMutsRegion(g, i, start, end)

		// Of all the possible silent mutations in the region, how many are
		// actually there?
		totalActual := 0
		for _, mut := range possibleMuts {
			if g.Nts[0][mut.Pos] == mut.To {
				totalActual++
			}
		}

		haveSite := hotspots.IsReSite(g.Nts[i][pos : pos+SITE_LEN])

		// If there aren't any actual mutations, the the probability that you
		// have a site here (if you got ths bit from this relative) is just 1
		// or 0.
		if totalActual == 0 {
			if haveSite {
				l.P += 1.0
			}
			continue
		}

		// This gives us an estimate of the probability of each nt in this
		// region getting mutated. Of course you don't know whether the muts
		// happened ages ago and the recombination yesterday though.
		perNt := float64(totalActual) / float64(len(possibleMuts))

		// Now look at the silent alternatives in the region.
		var env genomes.Environment
		env.Init(g, pos, SITE_LEN, 0)
		alternatives := env.FindAlternatives(6, true)

		var selected, total float64
		for _, alt := range alternatives {
			p := math.Pow(perNt, float64(alt.NumMuts))
			if hotspots.IsReSite(alt.Nts) != haveSite {
				selected += p
			}
			total += p
		}
		if haveSite {
			// If we have a site then we're talking about the probability of
			// losing it...
			l.P += 1.0 - (selected / total)
		} else {
			// ... and of gaining one if we don't
			l.P += selected / total
		}
	}
	// We're averaging these over all the relatives in the set.
	l.P /= float64(g.NumGenomes() - 1)
}

/*
Find anywhere either SC2 or one of the others has a site. Each Locus contains
an estimated probabilty that SC2 would have a site there after recombination
and point mutation from the family members in g
*/
func FindLoci(g *genomes.Genomes, window int) []Locus {
	ret := make([]Locus, 0)

	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < g.Length()-SITE_LEN; j++ {
			if hotspots.IsReSite(g.Nts[i][j : j+SITE_LEN]) {
				var locus Locus
				locus.Init(g, j, window)
				ret = append(ret, locus)
			}
		}
	}
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/Hassanin.fasta",
		"../fasta/WH1.orfs", false)
	g.RemoveGaps()
	loci := FindLoci(g, 50)
	for _, locus := range loci {
		fmt.Println(locus)
	}
}
