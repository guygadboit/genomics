package main

import (
	"fmt"
	"log"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"math/rand"
	"strings"
)

// What are the odds of finding a single random mut that matches the genomes in
// g (from 1 thru the end), but that doesn't match any in mask? Return the
// number of hits over its iterations. If silent, only look for silent
// mutations.
func OutgroupMontecarlo(g *genomes.Genomes, mask *genomes.Genomes,
	nd *mutations.NucDistro, its int, silent bool) int {
	hits := 0

	type mutation struct {
		pos   int
		newNt byte
	}
	seen := make(map[mutation]bool)

	for i := 0; i < its; i++ {
		var mut mutation
		for {
			mut.pos = rand.Intn(g.Length())
			existing := g.Nts[0][mut.pos]
			for {
				mut.newNt = nd.Random()
				if mut.newNt != existing {
					break
				}
			}

			if seen[mut] {
				continue
			}
			seen[mut] = true

			if silent {
				isSilent, _, err := genomes.IsSilentWithReplacement(g,
					mut.pos, 0, 0, []byte{mut.newNt})
				if err != nil && err.Error() == "Not in ORF" {
					// We are calling these silent (they don't change the
					// protein) in OutgroupMatches etc.
					isSilent = true
				}
				if !isSilent {
					continue
				}
			}
			break
		}

		// Now see if we match something in g, but not in mask
	genomes:
		for k := 1; k < g.NumGenomes(); k++ {
			if g.Nts[k][mut.pos] == mut.newNt {
				if mask != nil {
					for m := 0; m < mask.NumGenomes(); m++ {
						if mask.Nts[m][mut.pos] == mut.newNt {
							continue genomes
						}
					}
				}
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
	silent, nonSilent := 0, 0
	for i, match := range m {
		var tag string
		if showTag {
			tag = match.tag
		}
		var silence string
		if match.Silence == database.NON_SILENT {
			nonSilent++
		} else {
			silence = "*"
			silent++
		}
		s[i] = fmt.Sprintf("%c%d%c%s%s", match.From,
			match.Pos, match.To, silence, tag)
	}
	return strings.Join(s, ",") + fmt.Sprintf(" S:NS=%d:%d", silent, nonSilent)
}

func (m Matches) Contains(mut database.Mutation) bool {
	for _, match := range m {
		if match.Mutation == mut {
			return true
		}
	}
	return false
}

// If silent consider only muts that are SILENT or NOT_IN_ORF
func OutgroupMatches(g *genomes.Genomes,
	muts database.Mutations, tag string, show bool, silent bool) Matches {
	ret := make(Matches, 0)
	for _, m := range muts {
		if silent {
			switch m.Silence {
			case database.NON_SILENT:
				fallthrough
			case database.NOT_IN_ORF:
				fallthrough
			case database.UNKNOWN:
				continue
			}
		}
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
	OutgroupMatches(g, muts, "", true, false)
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

type Odds struct {
	Matches    int
	NonMatches int
}

type MatchOdds struct {
	Silent Odds
	All    Odds
}

/*
What are the odds of a mutation in genome 0 in g matching one of the others in
g, while not matching any in mask? Compute odds independently for silent and
non-silent matches
*/
func OutgroupOdds(g *genomes.Genomes, mask *genomes.Genomes) MatchOdds {
	var ret MatchOdds

	matchAny := func(nt byte, where *genomes.Genomes, from int, pos int) bool {
		for i := from; i < where.NumGenomes(); i++ {
			if where.Nts[i][pos] == nt {
				return true
			}
		}
		return false
	}

	for pos := 0; pos < g.Length(); pos++ {
		// Only consider mutations in coding sequences, as it's not clear how
		// often mutations outside coding sequences break things, or therefore
		// what their probabilities are. So we just skip them.
		_, _, err := g.Orfs.GetCodonOffset(pos)
		if err != nil {
			continue
		}

		// Now try all three possible mutations in this position
		for _, nt := range []byte{'G', 'C', 'A', 'T'} {
			if nt == g.Nts[0][pos] {
				// Not a mutation, so doesn't count as a match or non-match.
				// Just skip it.
				continue
			}

			// Does it match anything in g[1:] (the outgroup we're interested
			// in), while also not matching anything in the mask (the close
			// relatives we're masking out)?
			good := matchAny(nt, g, 1, pos) && !matchAny(nt, mask, 0, pos)

			if good {
				ret.All.Matches++
			} else {
				ret.All.NonMatches++
			}

			isSilent, _, err := genomes.IsSilentWithReplacement(g,
				pos, 0, 0, []byte{nt})
			if err != nil {
				log.Fatal(err)
			}

			if isSilent {
				if good {
					ret.Silent.Matches++
				} else {
					ret.Silent.NonMatches++
				}
			}
		}
	}
	return ret
}
