package main

import (
	"encoding/gob"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"genomics/outgroup"
	"genomics/utils"
	"math"
	"math/rand"
	"os"
	"path"
	"reflect"
	"slices"
)

const SITE_LEN = 6

// Represents anywhere you might end up with a site because somebody in the set
// has one.
type Locus struct {
	Pos       int
	SC2       bool    // Does SC2 have a site here?
	Count     int     // How many relatives have a site here
	SameCount int     // How many relatives have the same thing as SC2 here
	AltCount  int     // How many alternatives here
	P         float64 // Estimated probability of SC2 having a site here
}

func (l Locus) String() string {
	return fmt.Sprintf("%d %t %d %d %d %.04f",
		l.Pos+1, l.SC2, l.Count, l.SameCount, l.AltCount, l.P)
}

// How many unique silent mutations from og are possible?
func totalMuts(og []byte, alternatives genomes.Alternatives) int {
	var ret int
	for i, nt := range og {
		seen := make(map[byte]bool)
		for _, alt := range alternatives {
			altNt := alt.Nts[i]
			if seen[altNt] {
				continue
			}
			if altNt != nt {
				seen[altNt] = true
				ret++
			}
		}
	}
	return ret
}

// How many unique mutable positions in a set of mutations? (more than one
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

func (l *Locus) Init(g *genomes.Genomes,
	alternatives genomes.Alternatives,
	have0 bool, pos, window int,
	weightCb func(i, pos int, weight float64)) {

	// Put what SC2 has in as the first alternative so it's just a complete
	// set.
	newAlts := make(genomes.Alternatives, len(alternatives)+1)
	newAlts[0] = genomes.Alternative{0, g.Nts[0][pos : pos+SITE_LEN], 0}
	copy(newAlts[1:], alternatives)
	alternatives = newAlts

	l.Pos = pos
	l.SC2 = have0
	l.AltCount = len(alternatives)

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
		have := g.Nts[i][pos : pos+SITE_LEN]

		start := max(0, pos-window)
		end := min(g.Length(), pos+SITE_LEN+window)

		possBefore := mutations.PossibleSilentMutsRegion(g, i, start, pos)
		possAfter := mutations.PossibleSilentMutsRegion(g, i, pos+SITE_LEN, end)
		possibleMuts := make([]mutations.Mutation,
			0, len(possBefore)+len(possAfter))
		possibleMuts = append(possibleMuts, possBefore...)
		possibleMuts = append(possibleMuts, possAfter...)

		// Of all the possible silent mutations in the region, how many are
		// actually there?
		totalActual := 0
		for _, mut := range possibleMuts {
			if g.Nts[0][mut.Pos] == mut.To {
				totalActual++
			}
		}

		haveSite := hotspots.IsReSite(have)
		if haveSite {
			l.Count++
		}

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

		// pMut is an estimate of the probability of getting each possible
		// mutation.
		pMut := float64(totalActual) / float64(len(possibleMuts))

		// Downweight this result by how distant the relative is in this region
		weight := math.Pow(pMut, float64(totalActual))
		if weightCb != nil {
			weightCb(i, pos, weight)
		}

		// Now factor in the probability of keeping or losing the site due to
		// mutational drift.
		var p float64
		if haveSite {
			// If you have the site, then how many possible silent mutations
			// are there that would lose it. You need to miss them all.
			n := totalMuts(g.Nts[i][pos:pos+SITE_LEN], alternatives)
			p = math.Pow(1-pMut, float64(n))
		} else {
			// if you don't have the site, then how many mutations do you need
			// to get one? There might be more than one site in the
			// alternatives so we sum over them (quite often GAGACG and GAGACC
			// for example)
			if alternatives.Contains(have) {
				for _, alt := range alternatives {
					if hotspots.IsReSite(alt.Nts) {
						n := outgroup.CompareRegion(alt.Nts, have, 0, SITE_LEN)
						m := SITE_LEN - n

						// To get this site here, you need n changes, and m to
						// stay the same.
						p += math.Pow(pMut, float64(n)) *
							math.Pow(1-pMut, float64(m))
					}
				}
			}
		}
		l.P += p * weight
		totalWeight += weight
	}
	l.P /= totalWeight
}

type WeightMap struct {
	entries map[int]float64
	total   float64
}

func (wm *WeightMap) Init() {
	wm.entries = make(map[int]float64)
}

func (wm *WeightMap) Add(relative int, weight float64) {
	wm.entries[relative] += weight
	wm.total += weight
}

func (wm WeightMap) Summarize(g *genomes.Genomes) {
	type record struct {
		index  int
		weight float64
	}
	data := make([]record, 0, len(wm.entries))
	for k, v := range wm.entries {
		data = append(data, record{k, v})
	}
	utils.SortByKey(data, false, func(datum record) float64 {
		return datum.weight
	})
	for i := 0; i < min(20, len(data)); i++ {
		datum := data[i]
		fmt.Printf("%-20s %f\n", g.Names[datum.index], datum.weight/wm.total)
	}
}

type Weights struct {
	overall WeightMap
	perSite map[int]WeightMap
}

func (w *Weights) Init() {
	w.overall.Init()
	w.perSite = make(map[int]WeightMap)
}

func (w *Weights) Update(relative, pos int, weight float64) {
	var (
		perSite WeightMap
		there   bool
	)
	if perSite, there = w.perSite[pos]; !there {
		perSite.Init()
	}

	perSite.Add(relative, weight)
	w.perSite[pos] = perSite
	w.overall.Add(relative, weight)
}

/*
Find anywhere there is something that either is or could be silently mutated
into a site. Each Locus contains an estimated probability that SC2 would have a
site there after recombination and point mutation from the family members in g
*/
func FindLoci(g *genomes.Genomes, window int) ([]Locus, Weights) {
	ret := make([]Locus, 0)

	// Track which viruses contribute the most weight to the final result
	var weights Weights
	weights.Init()

	const (
		NO_SITE = iota
		HAVE
		HAVE_ALT
	)

	for i := 0; i < g.Length()-SITE_LEN; i++ {
		var status int

		if hotspots.IsReSite(g.Nts[0][i : i+SITE_LEN]) {
			status = HAVE
		}

		var env genomes.Environment
		err := env.Init(g, i, SITE_LEN, 0)
		if err != nil {
			// If we aren't in an ORF, then there are no silent alternatives
			// here.
			continue
		}
		alternatives := env.FindAlternatives(SITE_LEN, true)

		if status == NO_SITE {
			for _, alt := range alternatives {
				if hotspots.IsReSite(alt.Nts) {
					status = HAVE_ALT
					break
				}
			}
		}

		if status == NO_SITE {
			continue
		}

		var locus Locus
		locus.Init(g, alternatives, status == HAVE, i, window,
			func(i, pos int, weight float64) {
				weights.Update(i, pos, weight)
			})
		ret = append(ret, locus)
	}
	return ret, weights
}

func Montecarlo(g *genomes.Genomes,
	loci []Locus, iterations int) {
	r := rand.New(rand.NewSource(78))

	var totalSites int
	var good int
	for i := 0; i < iterations; i++ {
		var (
			maxGap   int
			numSites int
			lastPos  int
		)
		for i, locus := range loci {
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
			lastPos = locus.Pos
		}
		lastGap := g.Length() - lastPos
		maxGap = max(maxGap, lastGap)
		if numSites == 5 && maxGap < 8000 {
			good++
		}
		totalSites += numSites
	}

	passRate := float64(good) / float64(iterations)
	meanSites := float64(totalSites) / float64(iterations)
	fmt.Printf(`Montecarlo
==========
Average number of sites: %f
Fraction passing DEFUSE criteria: %f
`, meanSites, passRate)
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

func (w Weights) Summarize(g *genomes.Genomes) {
	fmt.Printf("Top contributing relatives by proportion of total weight:\n")
	fmt.Println("Overall")
	w.overall.Summarize(g)

	fmt.Println("Per site")
	sites := make([]int, 0, len(w.perSite))
	for k, _ := range w.perSite {
		sites = append(sites, k)
	}
	slices.Sort(sites)

	for _, pos := range sites {
		fmt.Printf("Location %d\n", pos+1)
		w.perSite[pos].Summarize(g)
	}
}

func main() {
	var (
		fname  string
		ref    int
		window int
		mc     bool
	)

	flag.IntVar(&ref, "ref", 0, "Use different genome for the reference")
	flag.StringVar(&fname, "fname", "HassaninPlus.fasta", "Alignment to use")
	flag.IntVar(&window, "window", 50, "What window size to use")
	flag.BoolVar(&mc, "mc", false, "Montecarlo")
	flag.Parse()

	g := genomes.LoadGenomes(path.Join("../fasta", fname),
		"../fasta/WH1.orfs", false)

	if ref != 0 {
		g = g.Swap(0, ref)
		g.RemoveGaps(true)
	}

	loci, weights := FindLoci(g, window)

	totalP := 0.0
	totalOff := 0.0
	actual := 0

	// Gotta do this before you sort the loci
	if mc {
		Montecarlo(g, loci, 100000)
	}

	utils.SortByKey(loci, true, func(l Locus) float64 {
		return l.P
	})

	correct := make([]float64, 0)

	for _, locus := range loci {
		fmt.Println(locus)
		if locus.SC2 {
			correct = append(correct, locus.P)
			actual++
		} else {
			totalOff += locus.P
		}
		totalP += locus.P
	}
	fmt.Printf(`
Number of sites expected: %f (actual %d)
Total of non-sites: %f

`, totalP, actual, totalOff)

	slices.Sort(correct)
	slices.Reverse(correct)

	fmt.Println("Prediction success of actual sites (per site and cumulative)")
	success := 1.0
	for i, p := range correct {
		success *= p
		fmt.Printf("%d: %f %f\n", i+1, p, success)
	}
	fmt.Println()

	weights.Summarize(g)
}
