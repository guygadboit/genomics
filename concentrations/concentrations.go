package main

import (
	"bufio"
	"encoding/gob"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/simulation"
	"genomics/stats"
	"log"
	"math"
	"math/rand"
	"os"
	"strings"
	"strconv"
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

type Concentrations struct {
	Genomes *genomes.Genomes
	Length  int
	MinMuts int
	Silent  bool
	Concs   []Concentration
}

func (c *Concentrations) Find(g *genomes.Genomes, length int, minMuts int,
	requireSilent bool) {
	c.Genomes = g
	c.Length = length
	c.MinMuts = minMuts
	c.Concs = findConcentrations(g, length, minMuts, requireSilent)
}

func (c *Concentrations) Save(fname string) {
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	enc := gob.NewEncoder(fp)
	err = enc.Encode(c)
	if err != nil {
		log.Fatal(err)
	}
}

func (c *Concentrations) Load(fname string) {
	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)
	dec := gob.NewDecoder(fp)
	err = dec.Decode(c)

	if err != nil {
		log.Fatal(err)
	}
}

/*
Subtract any that are at the same position (they're probably different
lengths-- you can use this to find the singles that aren't also in doubles)
*/
func (c *Concentrations) Subtract(other *Concentrations) {
	newConcs := make([]Concentration, 0)

	set := make(map[int]bool)
	for _, conc := range other.Concs {
		set[conc.Pos] = true
	}

	for _, conc := range c.Concs {
		if !set[conc.Pos] {
			newConcs = append(newConcs, conc)
		}
	}
	c.Concs = newConcs
}

func findConcentrations(g *genomes.Genomes, length int, minMuts int,
	requireSilent bool) []Concentration {
	ret := make([]Concentration, 0)

positions:
	for i := 0; i < g.Length()-length; i++ {
		for j := 0; j < g.NumGenomes(); j++ {
			for k := 0; k < j; k++ {
				_, silent, numMuts := genomes.IsSilent(g, i, length, j, k)
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

type GraphDatum struct {
	Key  string
	Real float64
	Sim  float64
}

func GraphData(fname string, realMap, simMap *TransitionMap) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatalf("Can't create file %s\n", fname)
	}
	defer f.Close()

	data := make(map[string]GraphDatum)

	l := NewCountList(realMap.Bidirection)
	SortCountList(l)

	for _, c := range l {
		k := c.Key.String()
		data[k] = GraphDatum{k,
			float64(realMap.Bidirection[c.Key]),
			float64(simMap.Bidirection[c.Key])}
	}

	var totalReal, totalSim float64
	for _, d := range data {
		totalReal += float64(d.Real)
		totalSim += float64(d.Sim)
	}

	for k, v := range data {
		data[k] = GraphDatum{k, v.Real / totalReal, v.Sim / totalSim}
	}

	fp := bufio.NewWriter(f)

	for _, c := range l {
		d := data[c.Key.String()]
		fmt.Fprintf(fp, "%s %.4f %.4f\n", d.Key, d.Real, d.Sim)
	}

	fp.Flush()
}

func CreateHighlights(concentrations []Concentration) []genomes.Highlight {
	ret := make([]genomes.Highlight, len(concentrations))
	for i, p := range concentrations {
		ret[i] = genomes.Highlight{p.Pos, p.Pos + p.Length, 'v'}
	}
	return ret
}

/*
Compare counts of concentrations to simulations. If iterations is -1, do an
exhaustive comparison. Otherwise do a Montecarlo with that many its. Returns
the real and simulated transition maps, and a contingency table of the average
numbers of concentrations in the real and simulated comparisons.
*/
func CompareToSim(g *genomes.Genomes, length int, minMuts int,
	requireSilent bool, iterations int,
	mutantFunc simulation.MutantFunc) (TransitionMap, TransitionMap,
	stats.ContingencyTable) {
	var simTotal, realTotal, silentTotal, numComparisons int
	n := g.NumGenomes()
	nd := mutations.NewNucDistro(g)
	var realMap, simMap TransitionMap

	realMap.Init()
	simMap.Init()

	comparePair := func(a, b int) {
		g2 := g.Filter(a, b)
		concs := findConcentrations(g2, length, minMuts, requireSilent)
		realMap.Combine(CountTransitions(g2, 0, 1, concs))
		realCount := len(concs)

		simG, numSilent := mutantFunc(g2, 0, 1, nd)
		concs = findConcentrations(simG, length, minMuts, requireSilent)
		simMap.Combine(CountTransitions(simG, 0, 1, concs))
		simCount := len(concs)

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

	var ct stats.ContingencyTable
	l := g.Length()

	average := func(count, total int) int {
		if total == 0 {
			return 0
		}
		return int(math.Round(float64(count) / float64(total)))
	}
	realAverage := average(realTotal, numComparisons)
	simAverage := average(simTotal, numComparisons)

	ct.Init(realAverage, l-realAverage, simAverage, l-simAverage)
	fmt.Println(ct)
	ct.FisherExact()

	/*
		fmt.Printf("%d.%d Average ratio: %.4f (%.2f silent muts)\n",
			length, minMuts, float64(realTotal)/float64(simTotal),
			float64(silentTotal)/float64(numComparisons))
	*/

	return realMap, simMap, ct
}

func main() {
	var requireSilent bool
	var iterations int
	var graphName string
	var simulateTriples bool
	var simulateTags bool
	var inputFile, orfs string
	var summary bool
	var outputClu string

	flag.BoolVar(&requireSilent, "silent", false, "Look at silent "+
		"(rather than all) mutations")
	flag.IntVar(&iterations, "its", 100, "Number of iterations")
	flag.StringVar(&graphName, "graph", "transitions.txt",
		"graph data filename")
	flag.StringVar(&inputFile, "input", "../fasta/SARS2-relatives.fasta",
		"Input file (aligned fasta)")
	flag.StringVar(&orfs, "orfs", "", "ORFS file")
	flag.BoolVar(&summary, "summary", false,
		"Just print a summary of the genomes")
	flag.StringVar(&outputClu, "clu", "", "Output a clu-style file of the "+
		"genomes specified. Use e.g. 0,1 for the first two")

	flag.BoolVar(&simulateTriples, "triples", false, "Look at triples")
	flag.BoolVar(&simulateTags, "tags", false, "Look at 6.4 tags")

	flag.Parse()

	if orfs == "" {
		log.Fatalf("This requires an ORFs file")
	}

	g := genomes.LoadGenomes(inputFile, orfs, false)
	g.RemoveGaps()

	if summary {
		mutations.Summary(g)
		return
	}

	if outputClu != "" {
		fields := strings.Split(outputClu, ",")
		which := make([]int, len(fields))
		for i, f := range fields {
			var err error
			which[i], err = strconv.Atoi(f)
			if err != nil {
				log.Fatalf("Invalid index specification: <%s>\n", outputClu)
			}
		}
		g2 := g.Filter(which...)
		var concs Concentrations
		concs.Find(g2, 2, 2, requireSilent)

		highlights := CreateHighlights(concs.Concs)
		g2.SaveWithTranslation("highlights.clu", highlights, which...)
		fmt.Printf("Written highlights.clu\n")
		return
	}


	/*
		var concs Concentrations
		//concs.Find(g, 2, 2, true)
		concs.Load("SARS1-Concs.gob")
	*/

	var f1, f2 simulation.MutantFunc
	if requireSilent {
		f1, f2 = simulation.MakeSimulatedMutant,
			simulation.MakeSimulatedMutant3
	} else {
		f1, f2 = simulation.MakeSimulatedMutant2,
			simulation.MakeSimulatedMutant4
	}

	realMap, simMap, ct := CompareToSim(g, 2, 2, requireSilent, iterations, f1)

	fmt.Println("Real transition map")
	realMap.Print()
	fmt.Println()

	fmt.Println("Sim transition map")
	simMap.Print()

	GraphData(graphName, &realMap, &simMap)
	fmt.Printf("Graph data written to %s\n", graphName)

	fmt.Printf("Frequency of doubles: OR=%.4f p=%g\n", ct.OR, ct.P)

	if simulateTriples {
		_, _, ct2 := CompareToSim(g, 3, 3, requireSilent, iterations, f2)
		fmt.Printf("Frequency of triples if you simulate "+
			"doubles: OR=%.4f p=%g\n", ct2.OR, ct2.P)
	}

	if simulateTags {
		_, _, ct3 := CompareToSim(g, 6, 4, requireSilent, iterations, f2)
		fmt.Printf("Frequency of 6.4 tags if you simulate "+
			"doubles: OR=%.4f p=%g\n", ct3.OR, ct3.P)
	}

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

	/*
		realMap, simMap := CompareToSim(g, 2, 2, true, 100, f)
		fmt.Println("Real transition map")
		realMap.Print()
		fmt.Println()

		fmt.Println("Sim transition map")
		simMap.Print()

		GraphData("transitions.txt", &realMap, &simMap)
	*/
	/*
		for _, c := range concs.Concs {
			fmt.Println(c)
		}
	*/

	/*
		var count, invertedCount int
		for i := 0; i < g.NumGenomes(); i++ {
			for j := 0; j < i; j++ {
				var singles, doubles Concentrations
				g2 := g.Filter(i, j)
				doubles.Find(g2, 2, 2, false)
				fmt.Printf("%d-%d: ", i, j)
				x, y := ShowDirections(g2, Transition{"CT", "TC"}, doubles.Concs)

				singles.Find(g2, 1, 1, false)
				singles.Subtract(&doubles)
				z, w := ShowDirections(g2, Transition{"C", "T"}, singles.Concs)

				inverted := (x < y) != (z < w)
				if inverted {
					fmt.Println("Inverted")
					invertedCount++
				}
				count++
			}
		}
		fmt.Printf("%d/%d are inverted (%.2f)\n", invertedCount, count,
			float64(invertedCount)/float64(count))
	*/
	/*
		i := 7
		for j := 0; j < g.NumGenomes(); j++ {
			var singles, doubles Concentrations
			g2 := g.Filter(i, j)
			doubles.Find(g2, 2, 2, false)
			fmt.Printf("%d-%d: ", i, j)
			ShowDirections(g2, Transition{"CT", "TC"}, doubles.Concs)

			singles.Find(g2, 1, 1, false)
			singles.Subtract(&doubles)
			ShowDirections(g2, Transition{"C", "T"}, singles.Concs)
		}
	*/

	return

	// TODO: Putting them both on the same graph would be nice, and you can use
	// that for KS testing externally as well. So output 3 columns
	/*
		realMap.GraphData("transitions.txt")
		simMap.GraphData("sim.txt")
	*/

	return
}
