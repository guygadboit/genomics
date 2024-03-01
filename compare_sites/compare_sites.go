package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"math"
	"net"
	"reflect"
	"slices"
)

type StringSet map[string]bool

func ToSet(st []string) StringSet {
	ret := make(StringSet)
	for _, s := range st {
		ret[s] = true
	}
	return ret
}

type Pattern struct {
	name  string
	which int // which genome from the original set this is
	nts   string
	pos   int
	count int // how many relatives it is in
}

type PatternSet map[Pattern]bool

// Maps each genome to its patterns
type Results map[string]PatternSet

func (r Results) Print() {
	for k, v := range r {
		fmt.Printf("%s\n", k)
		/*
			slices.SortFunc(v, func(a, b Pattern) int {
				return a.count - b.count
			})
		*/
		for k, _ := range v {
			fmt.Printf("%s %d %d\n", k.nts, k.pos, k.count)
		}
	}
}

func (p *Pattern) NumMuts(other *Pattern) int {
	var ret int
	for i := 0; i < len(p.nts); i++ {
		if p.nts[i] != other.nts[i] {
			ret++
		}
	}
	return ret
}

func ByLocation(r Results) map[int][]Pattern {
	ret := make(map[int][]Pattern)
	for _, s := range r {
		for k, _ := range s {
			ret[k.pos] = append(ret[k.pos], k)
		}
	}
	return ret
}

func mapPattern(g *genomes.Genomes, start, length int, results Results) {
	// Maps pattern to how many genomes have that pattern
	m := make(map[string]int)

	// Maps pattern to the indices of the genomes with that pattern at this
	// location.
	indices := make(map[string][]int)

	for i := 0; i < g.NumGenomes(); i++ {
		pat := string(g.Nts[i][start : start+length])
		m[pat]++
		indices[pat] = append(indices[pat], i)
	}

	for k, v := range m {
		if !utils.IsRegularPattern([]byte(k)) {
			continue
		}
		for _, index := range indices[k] {
			name := g.Names[index]
			ps, there := results[name]
			if !there {
				ps = make(PatternSet)
			}
			ps[Pattern{name, index, k, start, v}] = true
			results[name] = ps
		}
	}
}

/*
Look for anywhere any of the patterns appears in an alignment, or a
sequence that could be silently mutated into one, and see if any of those
are unique across the set
*/
func FindPatterns(g *genomes.Genomes, patterns []string) Results {
	ret := make(Results)

	patSet := ToSet(patterns)

	for i := 0; i < g.Length(); i++ {
		for j := 0; j < g.NumGenomes(); j++ {
			var env genomes.Environment
			err := env.Init(g, i, 6, j)
			if err != nil {
				continue
			}

			// Do we actually have the pattern here?
			sp := string(g.Nts[j][i : i+6])
			interested := patSet[sp]

			// Or is one of the silent alternatives the pattern?
			if !interested {
				for _, alt := range env.FindAlternatives(6) {
					if patSet[string(alt.Nts)] {
						interested = true
						break
					}
				}
			}

			if !interested {
				continue
			}

			mapPattern(g, i, 6, ret)
			break
		}
	}
	return ret
}

func mostDifferent(patterns []Pattern) (*Pattern, *Pattern, int) {
	mostMuts := 0
	var bestP, bestQ *Pattern

	for i, p := range patterns {
		for j, q := range patterns {
			nm := p.NumMuts(&q)
			if nm > mostMuts {
				mostMuts = nm
				bestP, bestQ = &patterns[i], &patterns[j]
			}
		}
	}
	return bestP, bestQ, mostMuts
}

type ContingencyTable struct {
	A, B, C, D int
	OR, p      float64
}

type FisherResult struct {
	OR, p float64
}

var FisherInput chan *ContingencyTable
var FisherOutput chan FisherResult

func FisherClient() {
	conn, err := net.Dial("unix", "./fisher.sock")
	if err != nil {
		log.Fatalf("Did you start calc_or.py?")
	}
	defer conn.Close()

	buf := make([]byte, 256)
	var result FisherResult

	for {
		ct := <-FisherInput
		msg := fmt.Sprintf("%s\n", ct.ToString())
		conn.Write([]byte(msg))
		_, err := conn.Read(buf)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Sscanf(string(buf), "%g %g", &result.OR, &result.p)
		FisherOutput <- result
	}
}

// We can work out the OR ourselves
func (c *ContingencyTable) CalcOR() float64 {
	aF, bF, cF, dF := float64(c.A), float64(c.B), float64(c.C), float64(c.D)
	c.OR = (aF / bF) / (cF / dF)
	if math.IsNaN(c.OR) || math.IsInf(c.OR, 1) || math.IsInf(c.OR, -1) {
		c.OR = 0
	}
	return c.OR
}

// But if we want to know the p-value we'll get that from our Python server.
func (c *ContingencyTable) DoFisherExact() {
	FisherInput <- c
	result := <-FisherOutput
	c.OR, c.p = result.OR, result.p
}

func (ct *ContingencyTable) ToString() string {
	return fmt.Sprintf("CT[%d %d %d %d]", ct.A, ct.B, ct.C, ct.D)
}

func FindSites(g *genomes.Genomes, which int, patterns []string) []int {
	ret := make([]int, 0)
	var s genomes.Search
	for _, pat := range patterns {
		for s.Init(g, which, []byte(pat), 0.0); !s.End(); s.Next() {
			pos, _ := s.Get()
			ret = append(ret, pos)
		}
	}
	return ret
}

// Returns which of the sites pos is in, or -1 if it's not in any of them
func InSites(sites []int, pos int) int {
	for i, site := range sites {
		start, end := site, site+6
		// You could do a binary search here but probably not necessary
		if pos >= start && pos < end {
			return i
		}
	}
	return -1
}

// Is there a silent mut between genomes a and b at pos?
func IsSilent(g *genomes.Genomes, a, b int, pos int) bool {
	aNts := g.Nts[a]
	bNts := g.Nts[b]

	if aNts[pos] == bNts[pos] {
		return false // not a mut at all
	}

	var envA, envB genomes.Environment
	err := envA.Init(g, pos, 1, a)
	if err != nil {
		return false
	}

	err = envB.Init(g, pos, 1, a)
	if err != nil {
		return false
	}

	return reflect.DeepEqual(envA.Protein(), envB.Protein())
}

// Make the contingency table they used in their preprint. Also return the
// maximum number of silent muts we found in one site between the two. If calcP
// find the p value, not just the OR (the p-value is much slower)
func SilentInPatterns(g *genomes.Genomes,
	a, b int, patterns []string, calcP bool) (ContingencyTable, int) {
	aNts := g.Nts[a]
	bNts := g.Nts[b]
	var ret ContingencyTable

	// tracks how many silent muts per site
	perSite := make(map[int]int)

	sites := FindSites(g, a, patterns)
	sites = append(sites, FindSites(g, b, patterns)...)
	slices.Sort(sites)

	for i := 0; i < g.Length(); i++ {
		if !utils.IsRegularNt(aNts[i]) {
			continue
		}

		if !utils.IsRegularNt(bNts[i]) {
			continue
		}

		silent := IsSilent(g, a, b, i)
		site := InSites(sites, i)
		in := site != -1

		if in {
			if silent {
				ret.A++
				perSite[site]++
			} else {
				ret.B++
			}
		} else {
			if silent {
				ret.C++
			} else {
				ret.D++
			}
		}
	}

	if calcP {
		ret.DoFisherExact()
	} else {
		ret.CalcOR()
	}

	// Find the max per site
	maxPerSite := 0
	for _, v := range perSite {
		if v > maxPerSite {
			maxPerSite = v
		}
	}

	return ret, maxPerSite
}

// The result of testing a pair
type Pair struct {
	a, b       int // Which two genomes these are
	ss         float64
	ct         ContingencyTable
	maxPerSite int
}

func TestAllPairs(g *genomes.Genomes, patterns []string) []Pair {
	ret := make([]Pair, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < i; j++ {
			ss := g.SequenceSimilarity(i, j) * 100
			ct, mps := SilentInPatterns(g, i, j, patterns, false)
			ret = append(ret, Pair{i, j, ss, ct, mps})
			fmt.Printf("%s/%s: (%.2f%% ss) %s\n",
				g.Names[i], g.Names[j], ss, ct.ToString())
		}
	}
	return ret
}

// Return where WH1 and RaTG13 rank in all the pairs
func RankPairs(pairs []Pair) (int, int) {
	wh1, rat := -1, -1

	/*
		slices.SortFunc(pairs, func(a, b Pair) int {
			return int(b.ct.OR - a.ct.OR)
		})
	*/

	slices.SortFunc(pairs, func(a, b Pair) int {
		if b.ct.p < a.ct.p {
			return 1
		} else {
			return -1
		}
	})

	for i, pair := range pairs {
		fmt.Printf("%d: %d vs %d OR=%f p=%g %s mps=%d ss=%.2f%%\n", i,
			pair.a, pair.b, pair.ct.OR, pair.ct.p, pair.ct.ToString(),
			pair.maxPerSite, pair.ss)
		if wh1 == -1 {
			if pair.a == 0 || pair.b == 0 {
				wh1 = i
			}
		}
		if rat == -1 {
			if pair.a == 463 || pair.b == 463 {
				rat = i
			}
		}
	}

	return wh1, rat
}

func init() {
	FisherInput = make(chan *ContingencyTable)
	FisherOutput = make(chan FisherResult)
	go FisherClient()
}

func main() {
	g := genomes.LoadGenomes("../fasta/more_relatives.fasta",
		"../fasta/WH1.orfs", false)
	/*
		g := genomes.LoadGenomes("../fasta/relatives.fasta",
			"../fasta/WH1.orfs", false)
	*/
	g.RemoveGaps()

	interesting := []string{
		"GGTCTC", // BsaI
		"GAGACC", // BsaI
		"CGTCTC", // BsmBI
		"GAGACG", // BsmBI
	}

	/*
	for i := 0; i < len(interesting); i++ {
		fmt.Printf("Trying %s\n", interesting[i])
		pairs := TestAllPairs(g, interesting[i:i+1], true)
		a, b := RankPairs(pairs)
		fmt.Printf("Ranks are %d,%d for %s\n", a, b, interesting[i])
	}
	*/

	// Controls, that code for LR, but aren't BsaI or BsmBI sites or anything
	/*
		interesting := []string{
			"TTACGC",
			"GCGTAA",
			"CTACGA",
			"GCGTAG",
		}
	*/

	// g.PrintSummary()
	TestAllPairs(g, interesting)
	return

	results := FindPatterns(g, interesting)
	interestingSet := ToSet(interesting)
	byLocation := ByLocation(results)

	for k, v := range byLocation {
		p, q, md := mostDifferent(v)

		for _, pat := range v {
			fmt.Printf("%d: %s has %s\n", k, pat.name, pat.nts)
		}

		if md >= 4 && (interestingSet[p.nts] || interestingSet[q.nts]) {
			ss := g.SequenceSimilarity(p.which, q.which) * 100
			ct, _ := SilentInPatterns(g, p.which, q.which, interesting, false)
			fmt.Printf("%d %s vs %s: %s/%s %d muts %.2f%% ss ", k,
				p.name, q.name, p.nts, q.nts, md, ss)
			fmt.Printf("CT[%d %d %d %d]\n", ct.A, ct.B, ct.C, ct.D)
		}
	}
}
