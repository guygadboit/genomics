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
	nd *mutations.NucDistro, its int) (int, int) {
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
	return hits, its
}

type Match struct {
	database.Mutation
	pangolins bool
	bats      bool
}

type Matches []Match

func (m Matches) ToString() string {
	s := make([]string, len(m))
	for i, match := range m {
		var pangolins, bats string
		if match.pangolins {
			pangolins = "p"
		}
		if match.bats {
			bats = "b"
		}
		s[i] = fmt.Sprintf("%c%d%c/%s%s", match.From, match.Pos, match.To,
			pangolins, bats)
	}
	return strings.Join(s, ",")
}

func OutgroupMatches(g *genomes.Genomes,
	muts database.Mutations, numBats int, show bool) Matches {
	ret := make(Matches, 0)
	for _, m := range muts {
		bats, pangolins := false, false
		for i := 1; i < g.NumGenomes(); i++ {
			if g.Nts[i][m.Pos-1] == m.To {
				if show {
					fmt.Printf("%s shared by %s\n", m.ToString(), g.Names[i])
				}
				if i-1 < numBats {
					bats = true
				} else {
					pangolins = true
				}
			}
		}
		if bats || pangolins {
			ret = append(ret, Match{m, pangolins, bats})
		}
	}
	return ret
}

func ShowOutgroupMatches(g *genomes.Genomes, muts database.Mutations) {
	OutgroupMatches(g, muts, 0, true)
}
