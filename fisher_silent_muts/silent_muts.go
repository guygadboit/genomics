package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
	"log"
)

type Positions map[int]bool

type Mutation struct {
	mutations.Mutation
	In	bool	// is this mutation inside the sites?
}

// Set mut.In for all the muts in mutations. You're "In" a site if there's a
// site there either in the 0th genome or in the which'th one. Return the
// positions of those that are in and out.
func CountInSites(muts []Mutation,
	g *genomes.Genomes, which int, sites [][]byte) (Positions, Positions) {

	in := make(Positions)
	out := make(Positions)

	// Index the muts by their positions, and store them all in "out"
	// initially. We will remove them from there as we find them in sites
	mutPositions := make(map[int]int)
	for i, mut := range muts {
		mutPositions[mut.Pos] = i
		out[mut.Pos] = true
		muts[i].In = false
	}

	handleMatch := func(s *genomes.Search, site []byte) {
		pos, err := s.Get()
		if err != nil {
			log.Fatal(err)
		}
		for i := 0; i < len(site); i++ {
			j, there := mutPositions[pos+i]
			if there {
				muts[j].In = true
				in[pos+i] = true
				delete(out, pos+i)
			}
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
	return in, out
}

func makeHighlights(actualIn, actualOut,
	possibleIn, possibleOut Positions) []genomes.Highlight {
	ret := make([]genomes.Highlight, 0)

	for pos, _ := range actualIn {
		ret = append(ret, genomes.Highlight{pos, pos+1, 'Z'})
	}

	for pos, _ := range actualOut {
		ret = append(ret, genomes.Highlight{pos, pos+1, 'z'})
	}

	for pos, _ := range possibleIn {
		if !actualIn[pos] {
			ret = append(ret, genomes.Highlight{pos, pos+1, 'X'})
		}
	}

	/*
	for pos, _ := range possibleOut {
		if !actualOut[pos] {
			ret = append(ret, genomes.Highlight{pos, pos+1, 'x'})
		}
	}
	*/

	return ret
}

func MonteCarlo(g *genomes.Genomes, muts []Mutation, its int) {
	for i := 0; i < its; i++ {
		sites := make([][]byte, 4)

		for j := 0; j < 2; j++ {
			sites[j] = utils.RandomNts(6)
			sites[j+2] = utils.ReverseComplement(sites[j])
		}
		TestGenomes(g, muts, sites)
	}
}

func TestGenomes(g *genomes.Genomes, muts []Mutation, sites [][]byte) {
	for i := 1; i < g.NumGenomes(); i++ {
		actual := make([]Mutation, 0)
		for _, mut := range muts {
			// All of the places where this relative actually has a silent mut
			if g.Nts[i][mut.Pos] == mut.To {
				actual = append(actual, mut)
			}
		}

		CountInSites(actual, g, i, sites)
		CountInSites(muts, g, i, sites)

		// Now the actual counts, which are based on mutations, not on
		// positions.
		var actualIn, actualOut int
		for _, mut := range actual {
			if mut.In {
				actualIn++
			} else {
				actualOut++
			}
		}

		var possibleIn, possibleOut int
		for _, mut := range muts {
			if mut.In {
				possibleIn++
			} else {
				possibleOut++
			}
		}

		var ct stats.ContingencyTable
		ct.Init(actualIn, actualOut, possibleIn, possibleOut)
		OR, p := ct.FisherExact()
		fmt.Printf("%s: OR=%f p=%g\n", g.Names[i], OR, p)
	}
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

	muts := make([]Mutation, 0)
	for _, mut := range mutations.PossibleSilentMuts(g, 0) {
		muts = append(muts, Mutation{mut, false})
	}

	fmt.Println("With the actual sites")
	TestGenomes(g, muts, sites)
	fmt.Println("MonteCarlo")
	MonteCarlo(g, muts, 1000)
}
