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

func (m Matches) ToString(showTag bool) string {
	s := make([]string, len(m))
	for i, match := range m {
		var tag string
		if showTag {
			tag = match.tag
		}
		s[i] = fmt.Sprintf("%c%d%c%s", match.From,
			match.Pos, match.To, tag)
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
	muts database.Mutations, tag string, show bool) Matches {
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
			ret = append(ret, Match{m, tag})
		}
	}
	return ret
}

func ShowOutgroupMatches(g *genomes.Genomes, muts database.Mutations) {
	OutgroupMatches(g, muts, "", true)
}

/*
If not strict, consider mutations to be "the same" if they're just in the same
place, even if they mutate to different things
*/
func RemoveIntersection(matchesA, matchesB Matches,
	strict bool) (Matches, Matches) {
	excludeA, excludeB := make(map[int]bool), make(map[int]bool)

	areSameStrict := func(a, b database.Mutation) bool {
		return a == b
	}

	areSameLoose := func(a, b database.Mutation) bool {
		return a.Pos == b.Pos
	}

	var areSame func(a, b database.Mutation) bool
	if strict {
		areSame = areSameStrict
	} else {
		areSame = areSameLoose
	}

	for i, a := range matchesA {
		for j, b := range matchesB {
			if areSame(a.Mutation, b.Mutation) {
				excludeA[i] = true
				excludeB[j] = true
			}
		}
	}

	filter := func(matches Matches, exclude map[int]bool) Matches {
		ret := make(Matches, 0)
		for i, m := range matches {
			if !exclude[i] {
				ret = append(ret, m)
			}
		}
		return ret
	}

	return filter(matchesA, excludeA), filter(matchesB, excludeB)
}
