package main

import (
	//"bufio"
	//"log"
	//"os"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/simulation"
	"math/rand"
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

type Transition struct {
	ANts string
	BNts string
}

func (t Transition) String() string {
	return fmt.Sprintf("%s<->%s", t.ANts, t.BNts)
}

/*
Return a copy in a "Canonical Order", so that we can use this as a map key on
both nt strings regardless of which is A and which is B
*/
func (t Transition) CanonicalOrder() Transition {
	var a, b string
	if t.ANts < t.BNts {
		a, b = t.ANts, t.BNts
	} else {
		a, b = t.BNts, t.ANts
	}
	return Transition{a, b}
}

type TransitionMap struct {
	Forwards  map[string]int // Counts by ANts
	Backwards map[string]int // Counts by BNts

	// Counts by all nts involved, regardless of direction
	Bidirection map[Transition]int
}

func (tm *TransitionMap) Init() {
	tm.Forwards = make(map[string]int)
	tm.Backwards = make(map[string]int)
	tm.Bidirection = make(map[Transition]int)
}

func (tm *TransitionMap) Add(t Transition) {
	tm.Forwards[t.ANts]++
	tm.Backwards[t.BNts]++
	tm.Bidirection[t.CanonicalOrder()]++
}

/*
Count which pairs of nts were involved the most often in transitions. Return an
array sorted by most frequent first.
*/
func CountTransitions(g *genomes.Genomes,
	a, b int, concentrations []Concentration) TransitionMap {
	var ret TransitionMap
	ret.Init()

	for _, conc := range concentrations {
		aNts := string(g.Nts[a][conc.Pos : conc.Pos+conc.Length])
		bNts := string(g.Nts[b][conc.Pos : conc.Pos+conc.Length])
		ret.Add(Transition{aNts, bNts})
	}

	return ret
}

func (tm *TransitionMap) Combine(other TransitionMap) {
	for k, v := range other.Forwards {
		tm.Forwards[k] += v
	}
	for k, v := range other.Backwards {
		tm.Backwards[k] += v
	}
	for k, v := range other.Bidirection {
		tm.Bidirection[k] += v
	}
}

/*
func (tm TransitionMap) Sort() TransitionList {
	ret := make(TransitionList, 0)
	for k, v := range tm {
		ret = append(ret, TransitionCount{k, v})
	}

	slices.SortFunc(ret, func(a, b TransitionCount) int {
		return b.Count - a.Count
	})

	return ret
}
*/

func (tm TransitionMap) Print() {
	fmt.Println("Forwards")

	clf := NewCountList(tm.Forwards)
	SortPrintCountList(clf)
	fmt.Println()

	fmt.Println("Backwards")
	clb := NewCountList(tm.Backwards)
	SortPrintCountList(clb)
	fmt.Println()

	fmt.Println("Bidirection")
	cli := NewCountList(tm.Bidirection)
	SortPrintCountList(cli)
	fmt.Println()
}

/*
func (tm TransitionMap) Total() int {
	var total int
	for _, v := range tm {
		total += v
	}
	return total
}

func (tm TransitionMap) GraphData(fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatalf("Can't create %s", fname)
	}
	defer f.Close()

	fp := bufio.NewWriter(f)

	total := float64(tm.Total())
	for _, tc := range tm.Sort() {
		fmt.Fprintf(fp, "%s<->%s %f\n", tc.ANts, tc.BNts,
			float64(tc.Count)/total)
	}

	fp.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func CreateHighlights(concentrations []Concentration) []genomes.Highlight {
	ret := make([]genomes.Highlight, len(concentrations))
	for i, p := range concentrations {
		ret[i] = genomes.Highlight{p.Pos, p.Pos + p.Length, 'v'}
	}
	return ret
}
*/

type MutantFunc func(*genomes.Genomes,
	int, int, *mutations.NucDistro) (*genomes.Genomes, int)

/*
Compare counts of concentrations to simulations. If iterations is -1, do an
exhaustive comparison. Otherwise do a Montecarlo with that many its. Returns
the real and simulated transition maps.
*/
func CompareToSim(g *genomes.Genomes, length int, minMuts int,
	requireSilent bool, iterations int,
	mutantFunc MutantFunc) (TransitionMap, TransitionMap) {
	var simTotal, realTotal, silentTotal, numComparisons int
	n := g.NumGenomes()
	nd := mutations.NewNucDistro(g)
	var realMap, simMap TransitionMap

	realMap.Init()
	simMap.Init()

	comparePair := func(a, b int) {
		g2 := g.Filter(a, b)
		concs := FindConcentrations(g2, length, minMuts, requireSilent)
		realMap.Combine(CountTransitions(g2, 0, 1, concs))
		realCount := len(concs)

		simG, numSilent := mutantFunc(g2, 0, 1, nd)
		concs = FindConcentrations(simG, length, minMuts, requireSilent)
		simMap.Combine(CountTransitions(simG, 0, 1, concs))
		simCount := len(concs)

		fmt.Printf("%d.%d %d-%d: %d %d %.2f (%d)\n",
			length, minMuts, a, b, realCount, simCount,
			float64(realCount)/float64(simCount), numSilent)

		if simCount > 0 && realCount > 0 {
			simTotal += simCount
			realTotal += realCount
			silentTotal += numSilent
			numComparisons++
		}
	}

	if iterations == -1 {
		for i := 0; i < n; i++ {
			for j := 0; j < i; j++ {
				comparePair(i, j)
			}
		}
	} else {
		for i := 0; i < iterations; i++ {
			var a, b int
			for {
				a, b = rand.Intn(n), rand.Intn(n)
				if a != b {
					break
				}
			}
			comparePair(a, b)
		}
	}

	fmt.Printf("%d.%d Average ratio: %.4f (%.2f silent muts)\n",
		length, minMuts, float64(realTotal)/float64(simTotal),
		float64(silentTotal)/float64(numComparisons))

	return realMap, simMap
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)

	/*
		g := genomes.LoadGenomes("../fasta/SARS1-relatives.fasta",
			"../fasta/SARS1.orfs", false)
	*/
	f := simulation.MakeSimulatedMutant
	g = g.Filter(5, 33)

	/*
		g = g.Filter(5, 33)
		concs := FindConcentrations(g, 2, 2, true)
		fmt.Printf("%d doubles\n", len(concs))
		transitions := CountTransitions(g, 0, 1, concs)
		tl := transitions.Sort()
		for _, t := range tl {
			t.Print()
		}
		highlights := CreateHighlights(concs)
		g.SaveWithTranslation("output.clu", highlights, 0, 1)
	*/

	realMap, simMap := CompareToSim(g, 2, 2, true, -1, f)
	fmt.Println("Real transition map")
	realMap.Print()
	fmt.Println()

	fmt.Println("Sim transition map")
	simMap.Print()

	// TODO: Putting them both on the same graph would be nice, and you can use
	// that for KS testing externally as well. So output 3 columns
	/*
		realMap.GraphData("transitions.txt")
		simMap.GraphData("sim.txt")
	*/

	return

	CompareToSim(g, 3, 3, true, 100, f)
	CompareToSim(g, 6, 4, true, 100, f)
}
