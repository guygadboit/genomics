package main

import (
	"encoding/gob"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"math"
	"math/rand"
	"os"
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
	P         float64 // Estimated probability of SC2 having a site here
}

func (l Locus) String() string {
	return fmt.Sprintf("%d %t %d %d %d %.04f",
		l.Pos+1, l.SC2, l.Count, l.SameCount, l.AltCount, l.P)
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

func countMuts(a, b []byte) int {
	var ret int
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			ret++
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

	var totalWeight float64
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

		var env genomes.Environment
		err := env.Init(g, pos, SITE_LEN, 0)
		if err != nil {
			// This can happen if you swapped genomes around and end up with -
			// in the alignment. In that case just skip this site.
			continue
		}
		alternatives := env.FindAlternatives(SITE_LEN, true)
		l.AltCount = len(alternatives) + 1

		// If there aren't any actual mutations, the probability that you have
		// a site here (if you got this bit from this relative) is just 1 or 0.
		// The weight is also 1, because no mutations.
		if totalActual == 0 {
			if haveSite {
				l.P += 1.0
			}
			totalWeight += 1.0
			continue
		}

		// perNt is an estimate of the probability of each nt in this region
		// getting mutated. Of course you don't know whether the muts happened
		// ages ago and the recombination yesterday though.
		numPoss := countPositions(possibleMuts)
		perNt := float64(totalActual) / float64(numPoss)

		// Downweight this result by how distant the relative is in this region
		weight := math.Pow(perNt, float64(totalActual))

		// Now factor in the probability of keeping or losing the site due to
		// mutational drift.
		var p float64
		if haveSite {
			slots := totalSlots(g.Nts[0][pos:pos+SITE_LEN], alternatives)
			// A mutation in any of those slots will silently remove the site.
			// Th probability of missing each one of them is 1-perNt, and so
			// the probability of missing all of them is 1-perNt ** slots
			p = math.Pow(1-perNt, float64(slots))
		} else {
			for _, alt := range alternatives {
				if haveSite {
					if hotspots.IsReSite(alt.Nts) {
						// To get the site you need to get lucky for each of
						// the mutations you need.
						p = math.Pow(perNt, float64(alt.NumMuts))
						break
					}
				}
			}
		}
		l.P += p * weight
		totalWeight += weight
	}
	l.P /= totalWeight
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

type Result struct {
	Count         int
	WeightedCount float64
}

type Results struct {
	TotalWeight float64
	TotalCount  int
	PerCount    map[int]Result
}

func (r *Results) Show() {
	for i := 5; i <= 8; i++ {
		if counts, there := r.PerCount[i]; there {
			rate := float64(counts.Count) / float64(r.TotalCount)
			weightedRate := counts.WeightedCount / r.TotalWeight
			fmt.Printf("%d: %.06f %.06g\n", i, rate, weightedRate)
		}
	}
}

func Montecarlo(g *genomes.Genomes,
	loci []Locus, iterations int) Results {
	r := rand.New(rand.NewSource(78))
	var ret Results
	ret.PerCount = make(map[int]Result)

	var totalSites int
	for i := 0; i < iterations; i++ {
		var (
			maxGap   int
			numSites int
			lastPos  int
			weight   float64
		)
		weight = 0.0
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
		lastGap := g.Length() - lastPos
		maxGap = max(maxGap, lastGap)
		if maxGap < 8000 {
			counts := ret.PerCount[numSites]
			counts.Count++
			counts.WeightedCount += weight
			ret.PerCount[numSites] = counts
		}
		totalSites += numSites
		ret.TotalWeight += weight
	}
	ret.TotalCount = iterations

	fmt.Printf("%.02f sites on average\n",
		float64(totalSites)/float64(iterations))

	return ret
}

func SaveLoci(loci []Locus) {
	fd, _ := os.Create("loci.gob")
	defer fd.Close()

	encoder := gob.NewEncoder(fd)
	encoder.Encode(loci)
}

func LoadLoci() []Locus {
	fd, err := os.Open("loci.gob")
	if err != nil {
		return nil
	}
	decoder := gob.NewDecoder(fd)
	ret := make([]Locus, 0)
	decoder.Decode(&ret)
	return ret
}

func main() {
	var (
		fname    string
		useSaved bool
		ref      int
	)

	flag.BoolVar(&useSaved, "load", false, "Use saved loci")
	flag.IntVar(&ref, "ref", 0, "Use different genome for the reference")
	flag.StringVar(&fname, "fname", "HassaninPlus.fasta", "Alignment to use")
	flag.Parse()

	g := genomes.LoadGenomes(path.Join("../fasta", fname),
		"../fasta/WH1.orfs", false)

	if ref != 0 {
		g = g.Swap(0, ref)
	}

	var loci []Locus
	if useSaved {
		loci = LoadLoci()
	} else {
		loci = FindLoci(g, 50)
		SaveLoci(loci)
	}

	for _, locus := range loci {
		fmt.Println(locus)
	}
}
