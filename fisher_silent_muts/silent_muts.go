package main

import (
	"fmt"
	"log"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/utils"
)

type Positions map[int]bool

// Return how many muts in muts are inside and outside the sites, wherever
// those sites appear either in the 0th or the which'th genome
func CountInSites(g *genomes.Genomes, which int,
	positions Positions) (int, int) {
	sites := [][]byte{
		[]byte("GGTCTC"),
		[]byte("GAGACC"),
		[]byte("CGTCTC"),
		[]byte("GAGACG"),
	}
	var inCount int

	// Make a copy since we destroy this. This starts off as the set of all
	// positions that are outside the sites.
	out := make(Positions)
	for k, _ := range positions {
		out[k] = true
	}

	handleMatch := func(s *genomes.Search, site []byte) {
		pos, err := s.Get()
		if err != nil {
			log.Fatal(err)
		}
		for i := 0; i < len(site); i++ {
			if out[pos+i] {
				out[pos+i] = false
				inCount++
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
	return inCount, len(out)
}

func main() {
	g := genomes.LoadGenomes("../fasta/CloseRelatives.fasta",
		"../fasta/WH1.orfs", false)

	mutations := mutations.PossibleSilentMuts(g, 0)
	possible := make(Positions)

	for _, mut := range mutations {
		possible[mut.Pos] = true
	}

	for i := 1; i < g.NumGenomes(); i++ {
		actual := make(Positions)
		for _, mut := range mutations {
			if g.Nts[1][mut.Pos] == mut.To {
				actual[mut.Pos] = true
			}
		}
		nsIn, nsOut := CountInSites(g, i, utils.Difference(possible, actual))
		sIn, sOut := CountInSites(g, i, actual)
		fmt.Println(nsIn, nsOut, sIn, sOut)
	}
}
