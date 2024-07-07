package main

import (
	"fmt"
	"math/rand"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
	"os"
	"bufio"
	"log"
)

type Positions map[int]bool

type Mutation struct {
	mutations.Mutation
	In bool // is this mutation inside the sites?
}

// Find all the positions which are inside a site either in the 0th or the
// which'th genome
func FindPositions(g *genomes.Genomes, which int, sites [][]byte) Positions {
	ret := make(Positions)
	handleMatch := func(s *genomes.Search, site []byte) {
		pos, err := s.Get()
		if err != nil {
			log.Fatal(err)
		}
		for i := 0; i < len(site); i++ {
			ret[pos+i] = true
		}
	}

	for _, site := range sites {
		var s genomes.Search
		for s.Init(g, 0, site, 0.0); !s.End(); s.Next() {
			handleMatch(&s, site)
		}
		for s.Init(g, which, site, 0.0); !s.End(); s.Next() {
			handleMatch(&s, site)
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
func SetIn(muts []Mutation, positions Positions) {
	for i, mut := range muts {
		muts[i].In = positions[mut.Pos]
	}
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

// Pass in either sites or positions.
func TestGenomes(g *genomes.Genomes,
	possible []Mutation, sites [][]byte,
	positions Positions, verbose bool) []Result {
	ret := make([]Result, 0)
	for i := 1; i < g.NumGenomes(); i++ {
		if positions == nil {
			positions = FindPositions(g, i, sites)
		}
		actual := FindActual(g, i, possible)
		SetIn(actual, positions)
		SetIn(possible, positions)
		ct := FindCT(actual, possible)

		OR, p := ct.FisherExact(stats.GREATER)

		if verbose {
			fmt.Println(ct)
			fmt.Printf("%s: OR=%f p=%g\n", g.Names[i], OR, p)
		}
		ret = append(ret, Result{i, OR, p})
	}
	return ret
}

// Make a similar set of positions to if you were looking for sites but just
// any old where. We do cluster them in groups of 6
func RandomPositions(g *genomes.Genomes, which int) Positions {
	ret := make(Positions)

	for i := 0; i < 8; i++ {
		start := rand.Intn(g.Length())
		for j := 0; j < 6; j++ {
			ret[start+j] = true
		}
	}

	return ret
}

func MonteCarlo(g *genomes.Genomes,
	possible []Mutation, its int, useSites bool) []Result {
	ret := make([]Result, 0)
	for i := 0; i < its; i++ {
		var positions Positions
		if useSites {
			sites := make([][]byte, 4)
			for j := 0; j < 2; j++ {
				sites[j] = utils.RandomNts(6)
				// sites[j+2] = utils.ReverseComplement(sites[j])
				sites[j+2] = utils.RandomNts(6)
			}
			positions = FindPositions(g, 1, sites)
		} else {
			positions = RandomPositions(g, 1)
		}
		ret = append(ret, TestGenomes(g, possible, nil, positions, false)...)
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

func main() {
	sites := [][]byte{
		[]byte("GGTCTC"),
		[]byte("GAGACC"),
		[]byte("CGTCTC"),
		[]byte("GAGACG"),
	}

	g := genomes.LoadGenomes("../fasta/CloseRelatives.fasta",
		"../fasta/WH1.orfs", false)

	possible := make([]Mutation, 0)
	for _, mut := range mutations.PossibleSilentMuts(g, 0) {
		possible = append(possible, Mutation{mut, false})
	}

	TestGenomes(g, possible, sites, nil, true)

	g = g.Filter(0, 1)
	mc := MonteCarlo(g, possible, 5000, true)
	OutputResults(mc, 1)
}
