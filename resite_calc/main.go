package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"genomics/outgroup"
	"math"
	"math/rand"
	"path"
	"reflect"
)

const SITE_LEN = 6

// Represents anywhere you might end up with a site because somebody in the set
// has one.
type Locus struct {
	Pos       int
	SC2       bool    // Does SC2 have a site here?
	Count     int     // How many relatives have a site here
	SameCount int     // How many relatives have the same thing as SC2 here
	AltCount  int     // How many silent alternatives here.
	P         float64 // Probability of SC2 having a site here after mutations
	SS        float64 // Best sequence similarity in window around site
}

func (l Locus) String() string {
	return fmt.Sprintf("%d %t %d %d %d %.04f %.04f",
		l.Pos+1, l.SC2, l.Count, l.SameCount, l.AltCount, l.P, l.SS)
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

func bestSS(g *genomes.Genomes, start, end, pos int) float64 {
	var ret float64

	for i := 1; i < g.NumGenomes(); i++ {
		// Exclude the site itself
		before := g.SubSequenceSimilarity(0, i, start, pos, true)
		after := g.SubSequenceSimilarity(0, i, pos, end, true)

		mean := (before + after) / 2
		if i == 1 {
			ret = mean
		} else {
			ret = max(ret, mean)
		}
	}
	return ret
}

func (l *Locus) Init(g *genomes.Genomes, pos int, window int) {
	l.Pos = pos
	l.SC2 = hotspots.IsReSite(g.Nts[0][pos : pos+SITE_LEN])

	for i := 1; i < g.NumGenomes(); i++ {
		if hotspots.IsReSite(g.Nts[i][pos : pos+SITE_LEN]) {
			l.Count++
		}
		if reflect.DeepEqual(g.Nts[i][pos:pos+SITE_LEN],
			g.Nts[0][pos:pos+SITE_LEN]) {
			l.SameCount++
		}
	}

	closest := outgroup.FindNumClosest(g, 0, pos,
		SITE_LEN, window, 3, outgroup.KEEP_BEST)
	other := closest[0].Which
	window = closest[0].Window

	start := max(0, pos-window)
	end := min(g.Length(), pos+window)
	l.SS = bestSS(g, start, end, pos)
	possibleMuts := mutations.PossibleSilentMutsRegion(g, other, start, end)

	// Of all the possible silent mutations in the region, how many are
	// actually there?
	totalActual := 0
	for _, mut := range possibleMuts {
		if g.Nts[0][mut.Pos] == mut.To {
			totalActual++
		}
	}

	haveSite := hotspots.IsReSite(g.Nts[other][pos : pos+SITE_LEN])
	if haveSite {
		l.Count++
	}

	var env genomes.Environment
	err := env.Init(g, pos, SITE_LEN, 0)
	if err != nil {
		// This can happen if you removed 0 from the alignment. We only do that
		// when we're doing the Montecarlo, which doesn't use P anyway.
		return
	}
	alternatives := env.FindAlternatives(SITE_LEN, true)
	l.AltCount = len(alternatives) + 1

	// If there aren't any actual mutations, the probability that you have
	// a site here (if you got this bit from this relative) is just 1 or 0.
	if totalActual == 0 {
		if haveSite {
			l.P = 1.0
		}
		return
	}

	// perNt is an estimate of the probability of each nt in this region
	// getting mutated. Of course you don't know whether the muts happened
	// ages ago and the recombination yesterday though.
	numPoss := countPositions(possibleMuts)
	perNt := float64(totalActual) / float64(numPoss)

	if haveSite {
		slots := totalSlots(g.Nts[0][pos:pos+SITE_LEN], alternatives)
		// A mutation in any of those slots will silently remove the site.
		// Th probability of missing each one of them is 1-perNt, and so
		// the probability of missing all of them is 1-perNt ** slots
		l.P = math.Pow(1-perNt, float64(slots))
	} else {
		for _, alt := range alternatives {
			if haveSite {
				if hotspots.IsReSite(alt.Nts) {
					// To get the site you need to get lucky for each of
					// the mutations you need.
					l.P = math.Pow(perNt, float64(alt.NumMuts))
					break
				}
			}
		}
	}
}

/*
Find anywhere either SC2 or one of the others has a site. Each Locus contains
an estimated probability that SC2 would have a site there after recombination
and point mutation from the family members in g
*/
func FindLoci(g *genomes.Genomes, window int) []Locus {
	ret := make([]Locus, 0)

	for i := 0; i < g.Length()-SITE_LEN; i++ {
		for j := 0; j < g.NumGenomes(); j++ {
			if hotspots.IsReSite(g.Nts[j][i : i+SITE_LEN]) {
				var locus Locus
				locus.Init(g, i, window)
				ret = append(ret, locus)
				break
			}
		}
	}
	return ret
}

func Montecarlo(g *genomes.Genomes,
	loci []Locus, iterations int) float64 {
	r := rand.New(rand.NewSource(78))

	var good, totalSites int
	for i := 0; i < iterations; i++ {
		var (
			maxGap   int
			numSites int
			lastPos  int
		)
		for i, locus := range loci {
			if r.Int()%g.NumGenomes() > locus.Count {
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
			lastPos = locus.Pos
		}
		maxGap = max(g.Length() - lastPos)
		if numSites == 5 && maxGap < 8000 {
			good++
		}
		totalSites += numSites
	}

	fmt.Printf("%.02f sites on average\n",
		float64(totalSites)/float64(iterations))

	return float64(good) / float64(iterations)
}

func main() {
	var (
		mc    bool
		fname string
	)

	flag.BoolVar(&mc, "mc", false, "Montecarlo chimeras")
	flag.StringVar(&fname, "fname", "HassaninPlus.fasta", "Alignment to use")
	flag.Parse()

	g := genomes.LoadGenomes(path.Join("../fasta", fname),
		"../fasta/WH1.orfs", false)
	g.RemoveGaps()

	if mc {
		// Exclude SC2 as it's causally downstream
		g = g.FilterOut(0)
	}

	loci := FindLoci(g, 50)
	var total float64

	if mc {
		good := Montecarlo(g, loci, 10000)
		fmt.Printf("Fraction of chimeras that are 6,<8k: %f\n",
			good)
	} else {
		for _, locus := range loci {
			fmt.Println(locus)
			total += locus.P
		}
		fmt.Printf("Total P: %f\n", total)
	}
}
