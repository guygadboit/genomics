package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
	"log"
	"math/rand"
	"os"
)

type Positions map[int]bool

type SitePositions struct {
	Starts Positions // Places where sites start
	All    Positions // Places that include any part of a site
}

type Mutation struct {
	mutations.Mutation
	In bool // is this mutation inside the sites?
}

// Find all the positions which are inside a site either in the 0th or the
// which'th genome. Look for sites "added" in 0 if added else "removed" in 0.
func FindPositions(g *genomes.Genomes,
	which int, sites [][]byte, added bool) SitePositions {
	var ret SitePositions
	ret.Starts = make(Positions)
	ret.All = make(Positions)
	handleMatch := func(s *genomes.Search, site []byte) {
		pos, err := s.Get()
		if err != nil {
			log.Fatal(err)
		}
		ret.Starts[pos] = true
		for i := 0; i < len(site); i++ {
			ret.All[pos+i] = true
		}
	}

	for _, site := range sites {
		var s genomes.Search
		if added {
			for s.Init(g, 0, site, 0.0); !s.End(); s.Next() {
				handleMatch(&s, site)
			}
		} else {
			for s.Init(g, which, site, 0.0); !s.End(); s.Next() {
				handleMatch(&s, site)
			}
		}
	}

	return ret
}

// Find the actual mutations between 0 and which given the possible ones
func FindActual(g *genomes.Genomes, which int, possible []Mutation) []Mutation {
	ret := make([]Mutation, 0)
	for _, mut := range possible {
		if g.Nts[which][mut.Pos] == mut.To {
			ret = append(ret, mut)
		}
	}
	return ret
}

// Set the In field on muts based on whether they are in positions
func SetIn(muts []Mutation, positions *SitePositions) {
	for i, mut := range muts {
		muts[i].In = positions.All[mut.Pos]
	}
}

// Find the contingency table based on a per-site rather than a per-nt analysis
func FindSiteCT(actual []Mutation, possible []Mutation,
	positions *SitePositions, length int) stats.ContingencyTable {
	var a, b, c, d int

	actualPositions := make(Positions)
	for _, mut := range actual {
		actualPositions[mut.Pos] = true
	}

	possiblePositions := make(Positions)
	for _, mut := range possible {
		possiblePositions[mut.Pos] = true
	}

	for i := 0; i < length-6; i++ {
		startsSite := positions.Starts[i]
		var hasActual, hasPossible bool
		for j := i; j < i+6; j++ {
			if actualPositions[j] {
				hasActual = true
			}
			if possiblePositions[j] {
				hasPossible = true
			}
		}
		if startsSite {
			if hasActual {
				a++
			}
			if hasPossible {
				c++
			}
		} else {
			if hasActual {
				b++
			}
			if hasPossible {
				d++
			}
		}
	}

	var ret stats.ContingencyTable
	ret.Init(a, b, c, d)
	return ret
}

func FindCT(actual []Mutation, possible []Mutation) stats.ContingencyTable {
	var actualIn, actualOut Positions
	var possibleIn, possibleOut Positions

	actualIn = make(Positions)
	actualOut = make(Positions)
	for _, mut := range actual {
		if mut.In {
			actualIn[mut.Pos] = true
		} else {
			actualOut[mut.Pos] = true
		}
	}

	possibleIn = make(Positions)
	possibleOut = make(Positions)
	for _, mut := range possible {
		if mut.In {
			possibleIn[mut.Pos] = true
		} else {
			possibleOut[mut.Pos] = true
		}
	}

	var ret stats.ContingencyTable
	ret.Init(len(actualIn), len(actualOut), len(possibleIn), len(possibleOut))
	return ret
}

type Result struct {
	genome int
	OR     float64
	p      float64
}

// Pass in either sites or positions (and nil for the other one)
func TestGenomes(g *genomes.Genomes,
	possible []Mutation, added bool, perSite bool, sites [][]byte,
	positions *SitePositions, verbose bool) []Result {
	ret := make([]Result, 0)
	for i := 1; i < g.NumGenomes(); i++ {
		if positions == nil {
			pos := FindPositions(g, i, sites, added)
			positions = &pos
		}
		actual := FindActual(g, i, possible)
		SetIn(actual, positions)
		SetIn(possible, positions)

		var ct stats.ContingencyTable
		if perSite {
			ct = FindSiteCT(actual, possible, positions, g.Length())
		} else {
			ct = FindCT(actual, possible)
		}

		OR, p := ct.FisherExact(stats.GREATER)

		if verbose {
			if added {
				fmt.Println("Sites added")
			} else {
				fmt.Println("Sites removed")
			}
			fmt.Println(ct)
			fmt.Printf("%s: OR=%f p=%g\n", g.Names[i], OR, p)
		}
		ret = append(ret, Result{i, OR, p})

		if verbose {
			highlights := makeHighlights(positions.All, 'v')
			fname := fmt.Sprintf("%d.clu", i)
			g.SaveWithTranslation(fname, highlights, 0, i)
			fmt.Printf("Wrote %s\n", fname)
		}
	}

	return ret
}

// Make a similar set of positions to if you were looking for sites but just
// any old where. We do cluster them in groups of 6
func RandomPositions(g *genomes.Genomes, which int) SitePositions {
	var ret SitePositions
	ret.All = make(Positions)
	ret.Starts = make(Positions)

	for i := 0; i < 35; i++ {
		start := rand.Intn(g.Length() - 6)
		ret.Starts[start] = true
		for j := 0; j < 6; j++ {
			ret.All[start+j] = true
		}
	}

	return ret
}

// What is the sequence similarity inside the site positions?
func SequenceSimilarity(g *genomes.Genomes, positions Positions) float64 {
	var same float64
	for pos, _ := range positions {
		if g.Nts[0][pos] == g.Nts[1][pos] {
			same++
		}
	}
	return same / float64(len(positions))
}

func MonteCarlo(g *genomes.Genomes,
	possible []Mutation, its int, useSites bool) []Result {
	ret := make([]Result, 0)
	for i := 0; i < its; i++ {
		var positions SitePositions
		if useSites {
			sites := make([][]byte, 4)
			for j := 0; j < 2; j++ {
				sites[j] = utils.RandomNts(6)
				sites[j+2] = utils.ReverseComplement(sites[j])
				// sites[j+2] = utils.RandomNts(6)
			}
			positions = FindPositions(g, 1, sites, true)
		} else {
			positions = RandomPositions(g, 1)
		}
		ret = append(ret, TestGenomes(g, possible,
			true, false, nil, &positions, false)...)
	}
	return ret
}

func OutputResults(results []Result, which int) {
	f, err := os.Create("ORs")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	orW := bufio.NewWriter(f)

	f2, err := os.Create("ps")
	if err != nil {
		log.Fatal(err)
	}
	defer f2.Close()
	pW := bufio.NewWriter(f2)

	for _, res := range results {
		if res.genome == which {
			fmt.Fprintln(orW, res.OR)
			fmt.Fprintln(pW, res.p)
		}
	}

	orW.Flush()
	pW.Flush()

	fmt.Println("Wrote ORs and ps")
}

func MakeTestGenomes(g *genomes.Genomes) {
	nd := mutations.NewNucDistro(g)

	var orfs genomes.Orfs
	ret := genomes.NewGenomes(orfs, 2)

	ret.Nts[0] = make([]byte, g.Length())
	ret.Nts[1] = make([]byte, g.Length())

	for i := 0; i < g.Length(); i++ {
		ret.Nts[0][i] = nd.Random()
		if rand.Float64() < 0.96 {
			ret.Nts[1][i] = ret.Nts[0][i]
		} else {
			ret.Nts[1][i] = nd.Random()
		}
	}
	ret.Names[0] = "Test Random Genome"
	ret.Names[1] = "Test Random Relative"

	ret.SaveMulti("random.fasta")
}

func makeHighlights(positions Positions, char byte) []genomes.Highlight {
	ret := make([]genomes.Highlight, 0)

	for pos, _ := range positions {
		ret = append(ret, genomes.Highlight{pos, pos + 1, char})
	}

	return ret
}

func main() {
	sites := [][]byte{
		[]byte("GGTCTC"),
		[]byte("GAGACC"),
		[]byte("CGTCTC"),
		[]byte("GAGACG"),
	}

	var fasta string
	var orfs string
	var added bool
	var doMC bool
	var perSite bool

	flag.StringVar(&fasta, "fasta",
		"../fasta/CloseRelatives.fasta", "relatives")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "orfs")
	flag.BoolVar(&added, "added",
		true, "Look at added rather than removed sites")
	flag.BoolVar(&doMC, "montecarlo", false, "Do the MonteCarlo")
	flag.BoolVar(&perSite, "per-site", false, "Look per site rather than per nt")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	possible := make([]Mutation, 0)
	for _, mut := range mutations.PossibleSilentMuts(g, 0) {
		possible = append(possible, Mutation{mut, false})
	}

	TestGenomes(g, possible, added, perSite, sites, nil, true)

	if doMC {
		g = g.Filter(0, 1)
		mc := MonteCarlo(g, possible, 5000, true)
		OutputResults(mc, 1)
	}
}
