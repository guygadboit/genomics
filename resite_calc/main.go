package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"math"
	"math/rand"
)

const SITE_LEN = 6

// Represents anywhere you might end up with a site because somebody in the set
// has one.
type Locus struct {
	Pos   int
	Count int     // How many relatives have a site here
	P     float64 // Probability of SC2 having a site here after mutations
}

// How many unique places can you have a silent mutation in?
func totalSlots(og []byte, alternatives genomes.Alternatives) int {
	var ret int
	for i, nt := range og {
		for _, alt := range alternatives {
			if alt.Nts[i] != nt {
				ret++
				break
			}
		}
	}
	return ret
}

// How many unique mutatable positions in a set of mutations? (more than one
// mut may be possible in the same place)
func countPositions(muts []mutations.Mutation) int {
	locations := make(map[int]bool)
	for _, mut := range muts {
		locations[mut.Pos] = true
	}
	return len(locations)
}

func (l *Locus) Init(g *genomes.Genomes, pos int, window int) {
	l.Pos = pos

	// Now compute P: an estimate of the probability that SC2 will end up with
	// a site here.
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
		if haveSite {
			l.Count++
		}

		// If there aren't any actual mutations, the probability that you have
		// a site here (if you got this bit from this relative) is just 1 or 0.
		if totalActual == 0 {
			if haveSite {
				l.P += 1.0
			}
			continue
		}

		// perNt is an estimate of the probability of each nt in this region
		// getting mutated. Of course you don't know whether the muts happened
		// ages ago and the recombination yesterday though.
		numPoss := countPositions(possibleMuts)
		perNt := float64(totalActual) / float64(numPoss)

		// Now look at the silent alternatives in the region.
		var env genomes.Environment
		env.Init(g, pos, SITE_LEN, 0)
		alternatives := env.FindAlternatives(SITE_LEN, true)

		if haveSite {
			slots := totalSlots(g.Nts[0][pos:pos+SITE_LEN], alternatives)
			// A mutation in any of those slots will silently remove the site.
			// Th probability of missing each one of them is 1-perNt, and so
			// the probability of missing all of them is 1-perNt ** slots
			l.P += math.Pow(1-perNt, float64(slots))
		} else {
			for _, alt := range alternatives {
				if haveSite {
					if hotspots.IsReSite(alt.Nts) {
						// To get the site you need to get lucky for each of
						// the mutations you need.
						l.P += math.Pow(perNt, float64(alt.NumMuts))
						break
					}
				}
			}
		}
	}
	// We're averaging these over all the relatives in the set.
	l.P /= float64(g.NumGenomes() - 1)
	fmt.Printf("At end: %f\n", l.P)
}

/*
Find anywhere either SC2 or one of the others has a site. Each Locus contains
an estimated probability that SC2 would have a site there after recombination
and point mutation from the family members in g
*/
func FindLoci(g *genomes.Genomes, window int) []Locus {
	ret := make([]Locus, 0)

	for i := 0; i < g.Length()-SITE_LEN; i++ {
		//i = 10443
		for j := 0; j < g.NumGenomes(); j++ {
			if hotspots.IsReSite(g.Nts[j][i : i+SITE_LEN]) {
				var locus Locus
				locus.Init(g, i, window)
				ret = append(ret, locus)
				break
			}
		}
		//break
	}
	return ret
}

func Montecarlo(loci []Locus, iterations int) float64 {
	r := rand.New(rand.NewSource(78))

	var (
		totalGood  float64
		totalSites int
	)
	for i := 0; i < iterations; i++ {
		var (
			maxGap   int
			numSites int
		)
		for i, locus := range loci {
			if r.Int()%2 != 0 {
				if r.Float64() > locus.P {
					continue
				}
				numSites++
				var gap int
				if i > 0 {
					gap = locus.Pos - loci[i-1].Pos
				} else {
					gap = locus.Pos
				}
				maxGap = max(maxGap, gap)
			}
		}
		totalSites += numSites
		if numSites == 5 && maxGap <= 8000 {
			totalGood++
		}
	}
	fmt.Printf("Average number of sites: %.02f\n",
		float64(totalSites)/float64(iterations))
	return totalGood / float64(iterations)
}

func main() {
	g := genomes.LoadGenomes("../fasta/Temmam.fasta",
		"../fasta/WH1.orfs", false)
	g.RemoveGaps()
	loci := FindLoci(g, 50)
	for _, locus := range loci {
		fmt.Println(locus)
	}
	result := Montecarlo(loci, 1000000)
	fmt.Printf("Result: %f\n", result)
}
