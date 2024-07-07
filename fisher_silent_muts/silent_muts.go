package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"math/rand"
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

// Make a similar set of positions to if you were looking for sites but just
// any old where.
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

func FindCT(possible []Mutation, actual []Mutation) stats.ContingencyTable {
	var a, b, c, d int
	for _, mut := range actual {
		if mut.In {
			a++
		} else {
			b++
		}
	}

	for _, mut := range possible {
		if mut.In {
			c++
		} else {
			d++
		}
	}

	var ret stats.ContingencyTable
	ret.Init(a, b, c, d)
	return ret
}

func TestGenomes(g *genomes.Genomes, possible []Mutation, sites [][]byte) {
	for i := 1; i < g.NumGenomes(); i++ {
		positions := FindPositions(g, i, sites)
		actual := FindActual(g, i, possible)
		SetIn(possible, positions)
		SetIn(actual, positions)
		ct := FindCT(possible, actual)
		
		OR, p := ct.FisherExact(stats.GREATER)
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

	possible := make([]Mutation, 0)
	for _, mut := range mutations.PossibleSilentMuts(g, 0) {
		possible = append(possible, Mutation{mut, false})
	}

	TestGenomes(g, possible, sites)

}
