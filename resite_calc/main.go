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
	"log"
	"math"
	"math/rand"
	"os"
	"path"
	"reflect"
	"slices"
	"strings"
)

const SITE_LEN = 6

var controlSites [][]byte

func makeControlSites() {
	controlSites = make([][]byte, 4)
	a := utils.RandomNts(6)
	b := utils.RandomNts(6)

	controlSites[0] = a
	controlSites[1] = utils.ReverseComplement(a)
	controlSites[2] = b
	controlSites[3] = utils.ReverseComplement(b)

	/*
		for _, site := range controlSites {
			fmt.Println(string(site))
		}
	*/
}

func isControlSite(pattern []byte) bool {
	for _, site := range controlSites {
		if reflect.DeepEqual(pattern, site) {
			return true
		}
	}
	return false
}

func isReSite(pattern []byte) bool {
	if controlSites != nil {
		return isControlSite(pattern)
	}
	return hotspots.IsReSite(pattern)
}

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
	α, β float64, have0 bool, pos int, window Window, censor bool,
	weightCb func(i, pos int, weight, pMut float64, numMuts int)) {

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
		if isReSite(g.Nts[i][pos : pos+SITE_LEN]) {
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
		win := window.Get(pos)

		start := max(0, pos-win)
		end := min(g.Length(), pos+SITE_LEN+win)

		var possibleMuts []mutations.Mutation
		if censor {
			possBefore := mutations.PossibleSilentMutsRegion(g, i, start, pos)
			possAfter := mutations.PossibleSilentMutsRegion(g, i, pos+SITE_LEN, end)
			possibleMuts = make([]mutations.Mutation,
				0, len(possBefore)+len(possAfter))
			possibleMuts = append(possibleMuts, possBefore...)
			possibleMuts = append(possibleMuts, possAfter...)
		} else {
			possibleMuts = mutations.PossibleSilentMutsRegion(g, i, start, end)
		}

		// Of all the possible silent mutations in the region, how many are
		// actually there?
		totalActual := 0
		for _, mut := range possibleMuts {
			if g.Nts[0][mut.Pos] == mut.To {
				totalActual++
			}
		}

		haveSite := isReSite(have)
		if haveSite {
			l.Count++
		}

		// pMut is an estimate of the probability of getting each possible
		// mutation.
		pMut := (float64(totalActual) + α) /
			(float64(len(possibleMuts)) + α + β)

		// If there aren't any actual mutations, the probability that you have
		// a site here (if you got this bit from this relative) is just 1 or 0.
		// The weight is also 1, because no mutations.
		if pMut == 0.0 {
			if haveSite {
				l.P += 1.0
			}
			if weightCb != nil {
				weightCb(i, pos, 1.0, pMut, 0)
			}
			totalWeight += 1.0
			continue
		}

		// Downweight this result by how distant the relative is in this region
		weight := math.Pow(pMut, float64(totalActual))
		if weightCb != nil {
			weightCb(i, pos, weight, pMut, totalActual)
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
					if isReSite(alt.Nts) {
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

type MutRate struct {
	rate    float64
	numMuts int
}

// We have an overall one of these and one per site
type WeightMap struct {
	// Maps relatives to weight
	entries  map[int]float64
	mutRates map[int]MutRate
	total    float64
}

func (wm *WeightMap) Init() {
	wm.entries = make(map[int]float64)
	wm.mutRates = make(map[int]MutRate)
}

func NewWeightMap() *WeightMap {
	var ret WeightMap
	ret.Init()
	return &ret
}

func (wm *WeightMap) Add(relative int, weight, pMut float64, numMuts int) {
	wm.entries[relative] += weight
	wm.total += weight
	wm.mutRates[relative] = MutRate{pMut, numMuts}
}

func (wm WeightMap) Summarize(g *genomes.Genomes, pos int) {
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

		name := strings.Clone(g.Names[datum.index])
		if pos != -1 &&
			isReSite(g.Nts[datum.index][pos:pos+SITE_LEN]) {
			name += "*"
		}

		pMut := wm.mutRates[datum.index]
		if pMut.rate == -1 {
			fmt.Printf("%-20s %f\n", name, datum.weight/wm.total)
		} else {
			fmt.Printf("%-20s %f %d %f\n", name, datum.weight/wm.total,
				pMut.numMuts, pMut.rate)
		}
	}
}

type Weights struct {
	overall WeightMap
	perSite map[int]*WeightMap
	window  Window
}

func (w *Weights) Init(window Window) {
	w.overall.Init()
	w.perSite = make(map[int]*WeightMap)
	w.window = window
}

func (w *Weights) Update(relative, pos int, weight, pMut float64, numMuts int) {
	if _, there := w.perSite[pos]; !there {
		w.perSite[pos] = NewWeightMap()
	}
	w.perSite[pos].Add(relative, weight, pMut, numMuts)
	w.overall.Add(relative, weight, -1, 0)
}

type Window interface {
	Get(pos int) int
}

type FixedWindow int

func (f FixedWindow) Get(pos int) int {
	return int(f)
}

type PerPos map[int]int

type AdaptiveWindow struct {
	g       *genomes.Genomes
	initial int
	algo    outgroup.Algo

	// Maps postions to window sizes
	cache map[int]int
}

func (aw *AdaptiveWindow) Init(g *genomes.Genomes,
	initial int, algo outgroup.Algo) {
	aw.g = g
	aw.initial = initial
	aw.algo = algo
	aw.cache = make(map[int]int)
}

func (aw *AdaptiveWindow) Get(pos int) int {
	if ret, there := aw.cache[pos]; there {
		return ret
	}
	prox := outgroup.FindNumClosest(aw.g, 0,
		pos, SITE_LEN, aw.initial, 3, aw.algo)
	ret := prox[0].Window

	aw.cache[pos] = ret
	return ret
}

/*
Find anywhere there is something that either is or could be silently mutated
into a site. Each Locus contains an estimated probability that SC2 would have a
site there after recombination and point mutation from the family members in g
*/
func FindLoci(g *genomes.Genomes,
	window Window, censor bool, α, β float64) ([]Locus, Weights) {
	ret := make([]Locus, 0)

	// Track which viruses contribute the most weight to the final result
	var weights Weights
	weights.Init(window)

	const (
		NO_SITE = iota
		HAVE
		HAVE_ALT
	)

	for i := 0; i < g.Length()-SITE_LEN; i++ {
		var status int

		if isReSite(g.Nts[0][i : i+SITE_LEN]) {
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
				if isReSite(alt.Nts) {
					status = HAVE_ALT
					break
				}
			}
		}

		if status == NO_SITE {
			continue
		}

		var locus Locus
		locus.Init(g, alternatives,
			α, β, status == HAVE, i, window, censor,
			func(i, pos int, weight, pMut float64, numMuts int) {
				weights.Update(i, pos, weight, pMut, numMuts)
			})
		ret = append(ret, locus)
	}
	return ret, weights
}

func Montecarlo(g *genomes.Genomes,
	loci []Locus, iterations int, verbose bool) float64 {
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

	if verbose {
		fmt.Printf(`Montecarlo
==========
Average number of sites: %f
Fraction passing DEFUSE criteria: %f
`, meanSites, passRate)
	}

	return passRate
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

type Counts struct {
	correct  []float64 // P for correctly predicted sites
	totalP   float64
	totalOff float64 // total P in places that aren't actually sites
	totalOn  float64 // total P in places that are
	actual   int     // number of actual sites
	totalNew float64 // total number of sites not shared by any originals
	success  float64 // total P for predicted sites / actual number
}

func (c *Counts) Init(loci []Locus) {
	c.correct = make([]float64, 0)
	for _, locus := range loci {
		if locus.SC2 {
			c.correct = append(c.correct, locus.P)
			c.totalOn += locus.P
			c.actual++
		} else {
			c.totalOff += locus.P
		}
		c.totalP += locus.P
		if locus.Count == 0 {
			c.totalNew += locus.P
		}
	}
	slices.Sort(c.correct)
	slices.Reverse(c.correct)
	c.success = c.totalOn / float64(c.actual)
}

func (w Weights) Summarize(g *genomes.Genomes) {
	fmt.Printf("Top contributing relatives by proportion of total weight:\n")
	fmt.Println("Overall")
	w.overall.Summarize(g, -1)

	fmt.Println("Per site")
	sites := make([]int, 0, len(w.perSite))
	for k, _ := range w.perSite {
		sites = append(sites, k)
	}
	slices.Sort(sites)

	for _, pos := range sites {
		pattern := g.Nts[0][pos : pos+SITE_LEN]
		s := fmt.Sprintf("Location %d (window %d) %s",
			pos+1, w.window.Get(pos), string(pattern))
		if isReSite(pattern) {
			s += "*"
		}
		fmt.Println(s)
		w.perSite[pos].Summarize(g, pos)
	}
}

func parseBeta(s string) (float64, float64) {
	if s == "" {
		return 0.0, 0.0
	}

	fields := strings.Split(s, ",")
	if len(fields) != 2 {
		log.Fatal("Invalid beta params <%s>", s)
	}

	α := utils.Atof(fields[0])
	β := utils.Atof(fields[1])
	return α, β
}

// The sample variance for picking n things out a population of N and getting a
// rate of p
func estimateVariance(p, n, N float64) (σ2 float64) {
	σ2 = ((p * (1 - p)) / n) * ((N - n) / (N - 1))
	return
}

// Pick Beta parameters based on your prior mean and variance
func pickBetaParams(μ, σ2 float64) (α float64, β float64) {
	t := (μ*(1-μ))/σ2 - 1
	α = μ * t
	β = (1 - μ) * t
	return
}

func ScanBetasTime(g *genomes.Genomes, its int) {
	for T := 1000; T >= 0; T -= 50 {
		μ := float64(T) * 1e-4
		σ2 := estimateVariance(μ, 100, 800)
		α, β := pickBetaParams(μ, σ2)
		win := FixedWindow(50)
		loci, _ := FindLoci(g, &win, true, α, β)
		result := Montecarlo(g, loci, its, false)
		fmt.Printf("T=%dy μ=%f σ2=%f -> Beta(%f, %f): %f\n",
			T, μ, σ2, α, β, result)
	}
}

// Pick Beta params based on a simple weighting
func pickBetaParamsWeight(μ, weight float64) (α float64, β float64) {
	α = weight
	β = weight * ((1 / μ) - 1)
	return
}

func ScanBetasRatio(g *genomes.Genomes, ratio float64, its int) {
	for i := 0; i < 10; i++ {
		f := 1 << i
		w := 1 / float64(f)
		α, β := pickBetaParamsWeight(ratio, w)
		win := FixedWindow(50)
		loci, _ := FindLoci(g, &win, true, α, β)
		result := Montecarlo(g, loci, its, false)
		fmt.Printf("Beta(%f, %f) w=%f ratio=%f %f\n", α, β, w, ratio, result)
	}
}

func countAll(g *genomes.Genomes, window Window, α, β float64) {
	for i := 0; i < g.NumGenomes(); i++ {
		g = g.Swap(0, i)
		loci, _ := FindLoci(g, window, true, α, β)
		var c Counts
		c.Init(loci)
		fmt.Printf("%s: %f %f\n", g.Names[0], c.success, c.totalOff)
		g = g.Swap(i, 0)
	}
}

func randomSites(g *genomes.Genomes, window Window, α, β float64) {
	for i := 0; i < 100; i++ {
		makeControlSites()
		loci, _ := FindLoci(g, window, true, α, β)
		var c Counts
		c.Init(loci)
		fmt.Printf("%f %f\n", c.success, c.totalOff)
	}
}

func main() {
	var (
		fname      string
		ref        int
		win        int
		grow       bool
		ogGrowAlgo bool
		mc         bool
		beta       string
		its        int
		scan       bool
		scanR      float64
		countall   bool
		randomize  bool
		censor     bool
	)

	flag.IntVar(&ref, "ref", 0, "Use different genome for the reference")
	flag.StringVar(&fname, "fname", "HassaninPlus.fasta", "Alignment to use")
	flag.IntVar(&win, "window", 50, "What window size to use")
	flag.BoolVar(&grow, "grow", false, "Grow the window to find 3 winners")
	flag.BoolVar(&ogGrowAlgo, "og-algo", false, "Use OG grow algo")
	flag.BoolVar(&mc, "mc", false, "Montecarlo")
	flag.BoolVar(&scan, "scan", false, "Explore Beta function params")
	flag.Float64Var(&scanR,
		"scanr", 0.0, "Try different Beta for given pseudocount ratio")

	// MBW using 1,79
	flag.StringVar(&beta, "beta", "", "Beta Function coefficients")
	flag.IntVar(&its, "iterations", 100000, "Montecarlo iterations")
	flag.BoolVar(&countall, "countall", false, "Count all")
	flag.BoolVar(&randomize, "random", false, "Try random control sites")
	flag.BoolVar(&censor, "censor", true, "Mask out the sites themselves")
	flag.Parse()

	g := genomes.LoadGenomes(path.Join("../fasta", fname),
		"../fasta/WH1.orfs", false)

	if ref != 0 {
		g = g.Swap(0, ref)
		g.RemoveGaps(true)
		g.SaveMulti("swapped.fasta")
		g.Orfs.Save("swapped-orfs.fasta")
		fmt.Println("New alignment saved to swapped.fasta " +
			"and swapped-orfs.fasta")
	}

	if scanR != 0.0 {
		ScanBetasRatio(g, scanR, its)
		return
	}

	if scan {
		ScanBetasTime(g, its)
		return
	}

	var window Window
	if grow {
		var aw AdaptiveWindow
		algo := outgroup.KEEP_BEST
		if ogGrowAlgo {
			algo = outgroup.ORIGINAL
		}
		aw.Init(g, win, algo)
		window = &aw
	} else {
		window = FixedWindow(win)
	}

	α, β := parseBeta(beta)

	if countall {
		countAll(g, window, α, β)
		return
	}

	if randomize {
		randomSites(g, window, α, β)
		return
	}

	loci, weights := FindLoci(g, window, censor, α, β)

	// Gotta do this before you sort the loci
	if mc {
		Montecarlo(g, loci, its, true)
	}

	utils.SortByKey(loci, true, func(l Locus) float64 {
		return l.P
	})

	for _, locus := range loci {
		fmt.Println(locus)
	}

	var c Counts
	c.Init(loci)

	fmt.Printf(`
Number of sites expected: %f (actual %d)
Total of non-sites: %f
Total of new sites: %f
Overall prediction success: %f

`, c.totalP, c.actual, c.totalOff, c.totalNew, c.success)

	fmt.Println("Prediction success of actual sites (per site and cumulative)")
	success := 1.0
	for i, p := range c.correct {
		success *= p
		fmt.Printf("%d: %f %f\n", i+1, p, success)
	}
	fmt.Println()

	weights.Summarize(g)
}
