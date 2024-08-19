package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"log"
	"os"
	"bufio"
)

type PosDatum struct {
	Possible   int    // Number of possible silent muts here
	Actual     []bool // Actual silent muts with each genome
	StartsSite []bool // Positions that start a site in each genome
	InSite     []bool // Positions that are in a site in each gnome

	// These aren't really needed for any of the code but nice to output to
	// match things up.
	Nt         []byte // The nts here in each genome
	Aa         []byte // The aa here in each genome
}

// We store one of these for each position, whether there's anything
// interesting there or not.
type PosInfo []PosDatum

func (p *PosDatum) Init(numGenomes int) {
	p.Actual = make([]bool, numGenomes)
	p.StartsSite = make([]bool, numGenomes)
	p.InSite = make([]bool, numGenomes)
	p.Nt = make([]byte, numGenomes)
	p.Aa = make([]byte, numGenomes)
}

// map of positions to possible mutations at that position
type PossibleMap map[int][]mutations.Mutation

func NewPossibleMap(muts []mutations.Mutation) PossibleMap {
	ret := make(map[int][]mutations.Mutation)
	for _, mut := range muts {
		_, there := ret[mut.Pos]
		if !there {
			ret[mut.Pos] = []mutations.Mutation{mut}
		} else {
			ret[mut.Pos] = append(ret[mut.Pos], mut)
		}
	}
	return ret
}

/*
For all the genomes and all the sites, annotate each position with whether it
starts a site, is in a site, and how many possible and actual mutations there
are there with each genome.
*/
func FindPositionInfo(g *genomes.Genomes,
	possible PossibleMap, sites [][]byte) PosInfo {
	ret := make(PosInfo, g.Length())

	translations := make([]genomes.TranslationMap, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		translations[i] = genomes.NewTranslationMap(genomes.Translate(g, i))
	}

	for i := 0; i < g.Length(); i++ {
		muts := possible[i]

		ret[i].Init(g.NumGenomes())
		ret[i].Possible = len(muts)

		for j := 1; j < g.NumGenomes(); j++ {
			for _, mut := range muts {
				if g.Nts[j][mut.Pos] == mut.To {
					ret[i].Actual[j] = true
				}
			}
		}

		for j := 0; j < g.NumGenomes(); j++ {
			ret[i].Nt[j] = g.Nts[j][i]
			aa := translations[j][i].Aa
			if aa == 0 {
				aa = '-'
			}
			ret[i].Aa[j] = aa
		}
	}

	handleMatch := func(s *genomes.Search, which int, site []byte) {
		pos, err := s.Get()
		if err != nil {
			log.Fatal(err)
		}
		ret[pos].StartsSite[which] = true
		for i := 0; i < len(site); i++ {
			ret[pos+i].InSite[which] = true
		}
	}

	for i := 0; i < g.NumGenomes(); i++ {
		for _, site := range sites {
			var s genomes.Search
			for s.Init(g, i, site, 0.0); !s.End(); s.Next() {
				handleMatch(&s, i, site)
			}
		}
	}

	return ret
}

func (p *PosInfo) SaveTSV() {
	f, err := os.Create("posinfo.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	fmt.Fprintf(w, "pos\tposs\t")
	numGenomes := len((*p)[0].Actual)
	for i := 0; i < numGenomes; i++ {
		fmt.Fprintf(w, "nt-%d\taa-%d\tact-%d\tstart-%d\tin-%d\t", i, i, i, i, i)
	}
	fmt.Fprintf(w, "\n")

	for i, datum := range *p {
		fmt.Fprintf(w, "%d\t%d\t", i, datum.Possible)
		for j := 0; j < numGenomes; j++ {
			fmt.Fprintf(w, "%c\t", datum.Nt[j])
			fmt.Fprintf(w, "%c\t", datum.Aa[j])
			fmt.Fprintf(w, "%t\t", datum.Actual[j])
			fmt.Fprintf(w, "%t\t", datum.StartsSite[j])
			fmt.Fprintf(w, "%t\t", datum.InSite[j])
		}
		fmt.Fprintf(w, "\n")
	}
	w.Flush()
	fmt.Printf("Wrote posinfo.tsv\n")
}

func (p *PosInfo) ShowSites() {
	for i := 0; i < len((*p)[0].StartsSite); i++ {
		fmt.Printf("Genome %d\n", i)
		for pos, pd := range *p {
			if pd.StartsSite[i] {
				fmt.Printf("Site at %d (0-based)\n", pos)
			}
		}
	}
}
