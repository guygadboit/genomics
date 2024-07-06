package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"log"
)

type Positions map[int]bool

// Return how the positions of muts inside and outside the sites, wherever
// those sites appear either in the 0th or the which'th genome
func CountInSites(g *genomes.Genomes, which int,
	positions Positions) (Positions, Positions) {
	sites := [][]byte{
		[]byte("GGTCTC"),
		[]byte("GAGACC"),
		[]byte("CGTCTC"),
		[]byte("GAGACG"),
	}
	// Make a copy since we destroy this. This starts off as the set of all
	// positions that are outside the sites.
	out := make(Positions)
	for k, _ := range positions {
		out[k] = true
	}

	in := make(Positions)

	handleMatch := func(s *genomes.Search, site []byte) {
		pos, err := s.Get()
		if err != nil {
			log.Fatal(err)
		}
		for i := 0; i < len(site); i++ {
			if out[pos+i] {
				out[pos+i] = false
				in[pos+i] = true
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
		ret = append(ret, genomes.Highlight{pos, pos+1, 'L'})
	}

	for pos, _ := range actualOut {
		ret = append(ret, genomes.Highlight{pos, pos+1, 'l'})
	}

	for pos, _ := range possibleIn {
		ret = append(ret, genomes.Highlight{pos, pos+1, 'X'})
	}

	for pos, _ := range possibleOut {
		ret = append(ret, genomes.Highlight{pos, pos+1, 'x'})
	}

	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/CloseRelatives.fasta",
		"../fasta/WH1.orfs", false)

	mutations := mutations.PossibleSilentMuts(g, 0)

	// All of the places where it's possible to have a silent mut
	possible := make(Positions)

	for _, mut := range mutations {
		possible[mut.Pos] = true
	}

	for i := 1; i < g.NumGenomes(); i++ {
		actual := make(Positions)
		for _, mut := range mutations {
			// All of the places where this relative actually has a silent mut
			if g.Nts[i][mut.Pos] == mut.To {
				actual[mut.Pos] = true
			}
		}
		possibleIn, possibleOut := CountInSites(g, i, possible)
		actualIn, actualOut := CountInSites(g, i, actual)

		var ct stats.ContingencyTable
		ct.Init(len(actualIn), len(actualOut),
			len(possibleIn), len(possibleOut))
		OR, p := ct.FisherExact()
		fmt.Println(ct)
		fmt.Printf("%s: OR=%f p=%f\n", g.Names[i], OR, p)

		highlights := makeHighlights(actualIn, actualOut,
			possibleIn, possibleOut)
		g.SaveWithTranslation(fmt.Sprintf("%d.clu", i), highlights, 0, i)
	}
}
