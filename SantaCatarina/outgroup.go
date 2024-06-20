package main

import (
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"log"
	"math/rand"
	"fmt"
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

// Return which muts match genomes in g from "from" onwards
func OutgroupMatches(g *genomes.Genomes,
	from int, muts database.Mutations) database.Mutations {
	ret := make(database.Mutations, 0)
mutations:
	for _, m := range muts {
		for i := from; i < g.NumGenomes(); i++ {
			if g.Nts[i][m.Pos-1] == m.To {
				ret = append(ret, m)
				continue mutations
			}
		}
	}
	return ret
}

/*
If not strict, consider mutations to be "the same" if they're just in the same
place, even if they mutate to different things
*/
func RemoveIntersection(matchesA, matchesB database.Mutations,
	strict bool) (database.Mutations, database.Mutations) {
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
			if areSame(a, b) {
				excludeA[i] = true
				excludeB[j] = true
			}
		}
	}

	filter := func(matches database.Mutations,
		exclude map[int]bool) database.Mutations {
		ret := make(database.Mutations, 0)
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
				// We're already skipping anything outside an ORF, so don't
				// expect an error.
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

func ShowAllPossibleSilentMuts(g *genomes.Genomes) {
	for pos := 0; pos < g.Length(); pos++ {
		nt := g.Nts[0][pos]
		for _, replacement := range []byte{'G', 'A', 'C', 'T'} {
			if replacement == nt {
				continue
			}
			silent, _, err := genomes.IsSilentWithReplacement(g, pos, 0, 0,
				[]byte{replacement})
			if err != nil {
				continue
			}
			if silent {
				fmt.Printf("%c%d%c\n", nt, pos+1, replacement)
			}
		}
	}
}
