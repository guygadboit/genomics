package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"math/rand"
	"strings"
)

// What are the odds of finding a single random mut in the genomes in g (from
// 1 thru the end)?
func OutgroupMontecarlo(g *genomes.Genomes,
	nd *mutations.NucDistro, its int) int {
	hits := 0

	for i := 0; i < its; i++ {
		pos := rand.Intn(g.Length())

		existing := g.Nts[0][pos]
		var newNt byte
		for {
			newNt = nd.Random()
			if newNt != existing {
				break
			}
		}

		for k := 1; k < g.NumGenomes(); k++ {
			if g.Nts[k][pos] == newNt {
				hits++
				break
			}
		}
	}
	return hits
}

type Match struct {
	database.Mutation
	tag string
}

type Matches []Match

func (m Matches) ToString() string {
	s := make([]string, len(m))
	for i, match := range m {
		s[i] = fmt.Sprintf("%c%d%c%s", match.From,
			match.Pos, match.To, match.tag)
	}
	return strings.Join(s, ",")
}

func (m Matches) Contains(mut database.Mutation) bool {
	for _, match := range m {
		if match.Mutation == mut {
			return true
		}
	}
	return false
}

func OutgroupMatches(g *genomes.Genomes,
	muts database.Mutations, tag string, show bool, exclude Matches) Matches {
	ret := make(Matches, 0)
	for _, m := range muts {
		found := false
		for i := 1; i < g.NumGenomes(); i++ {
			if g.Nts[i][m.Pos-1] == m.To {
				if show {
					fmt.Printf("%s shared by %s\n", m.ToString(), g.Names[i])
				}
				found = true
			}
		}
		if found {
			if exclude != nil && exclude.Contains(m) {
				continue
			}
			ret = append(ret, Match{m, tag})
		}
	}
	return ret
}

func ShowOutgroupMatches(g *genomes.Genomes, muts database.Mutations) {
	OutgroupMatches(g, muts, "", true, nil)
}
