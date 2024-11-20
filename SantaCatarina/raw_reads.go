package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
	"genomics/utils"
	"slices"
	"strings"
)

/*
Represents a match between a subpopulation of one genome and WH1
*/
type Match struct {
	genome        int
	numMatches    int                 // the number of matches in the reads
	opportunities int                 // number of places where genome differs from WH1
	positions     []utils.OneBasedPos // where the matches are
}

func (m *Match) Ratio() float64 {
	return float64(m.numMatches) / float64(m.opportunities)
}

/*
Return a new pileup containing only the matches that don't match the reference
and that have a depth >= minDepth.
*/
func Alternatives(ref *genomes.Genomes,
	pu *pileup.Pileup, minDepth int) *pileup.Pileup {
	var ret pileup.Pileup
	ret.Init()

	for i, _ := range pu.Records {
		rec := &pu.Records[i]
		reads := make([]pileup.Read, 0)
		for _, read := range rec.Reads {
			if ref.Nts[0][rec.Pos] == read.Nt {
				continue
			}
			if read.Depth < minDepth {
				continue
			}
			reads = append(reads, read)
		}
		if len(reads) != 0 {
			ret.Add(rec.Pos, reads)
		}
	}
	return &ret
}

/*
Do any of the secondary matches in a bunch of pileups match each other? Just
print them all out for now. other is the index of the "other" genome you are
looking for whether there is a QS resembling it.
*/
func MatchMulti(ref *genomes.Genomes,
	minDepth int, other int, pileups ...*pileup.Pileup) {
	apus := make([]*pileup.Pileup, len(pileups))
	for i, p := range pileups {
		apus[i] = Alternatives(ref, p, minDepth)
	}

	for pos := 0; pos < ref.Length(); pos++ {
		var lineStarted bool
		for _, alt := range apus {
			rec := alt.Get(pos)
			if rec == nil {
				continue
			}
			refNt := ref.Nts[0][pos]
			otherNt := ref.Nts[other][pos]
			different := refNt != otherNt

			if !lineStarted {
				var highlight string
				if different {
					highlight = "*"
				}
				fmt.Printf("%5d %c%c%s: ", pos, refNt, otherNt, highlight)
				lineStarted = true
			}
			for _, read := range rec.Reads {
				var special string
				if read.Nt == otherNt {
					special = "+"
				}
				fmt.Printf("%c%s", read.Nt, special)
			}
			fmt.Printf("-")
		}
		if lineStarted {
			fmt.Printf("\n")
		}
	}
}

func makeString(p []utils.OneBasedPos) string {
	s := make([]string, len(p))
	for i, pos := range p {
		s[i] = fmt.Sprintf("%d", pos)
	}
	return strings.Join(s, " ")
}

func MatchReads(g *genomes.Genomes, pu *pileup.Pileup, minDepth int) {
	alt := Alternatives(g, pu, minDepth)

	counts := make(map[int]int)
	totals := make(map[int]int)
	positions := make(map[int][]utils.OneBasedPos)

	for i := 0; i < g.Length(); i++ {
		rec := alt.Get(i)
		if rec == nil {
			continue
		}
		for j := 1; j < g.NumGenomes(); j++ {
			if g.Nts[j][i] == g.Nts[0][i] {
				continue
			}
			for _, read := range rec.Reads {
				totals[j]++
				if read.Nt == g.Nts[j][i] {
					counts[j]++
					positions[j] = append(positions[j], utils.OneBasedPos(i+1))
				}
			}
		}
	}

	matches := make([]Match, g.NumGenomes()-1)
	for i := 1; i < g.NumGenomes(); i++ {
		matches[i-1] = Match{i, counts[i], totals[i], positions[i]}
	}
	slices.SortFunc(matches, func(a, b Match) int {
		if a.Ratio() < b.Ratio() {
			return 1
		}
		if a.Ratio() > b.Ratio() {
			return -1
		}
		return 0
	})

	for _, m := range matches {
		fmt.Printf("%d: %d/%d %s %.4f\n",
			m.genome, m.numMatches,
			m.opportunities, makeString(m.positions), m.Ratio())
	}
}
