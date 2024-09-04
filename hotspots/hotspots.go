package hotspots

import (
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"math/rand"
)

var RE_SITES = [][]byte{
	[]byte("GGTCTC"), // BsaI
	[]byte("GAGACC"), // BsaI
	[]byte("CGTCTC"), // BsmBI
	[]byte("GAGACG"), // BsmBI
}

type Where int

const (
	SITE_IN_A Where = iota
	SITE_IN_B
	SITE_IN_EITHER
)

type Mutation struct {
	mutations.Mutation
	In0 bool // is this mutation inside the sites in the first genome?
	In1 bool // is this mutation inside the sites in the second genome?
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

type CalcCT struct {
	/*
		Calculate a Contingency Table in a variety of correct or incorrect ways
		given the actual mutations, the possible mutations, where the sites are,
		and the length of the genome. If correctDoubles then count where sites
		appear in both genomes correctly (only applies when where ==
		SITE_IN_EITHER)
	*/
	Calc func(posInfo PosInfo, which int,
		where Where, correctDoubles bool) stats.ContingencyTable

	// How long a window to use when looking at actually and possibly silently
	// mutated sequences.
	Window int
}

// How many "hits" to count for a match.
func findScore(where Where, correctDoubles, inA, inB bool) int {
	var score int

	switch where {
	case SITE_IN_A:
		if inA {
			score = 1
		}
	case SITE_IN_B:
		if inB {
			score = 1
		}
	case SITE_IN_EITHER:
		if inA || inB {
			score = 1
		}
		if correctDoubles && inA && inB {
			score = 2
		}
	}
	return score
}

func findCT(posInfo PosInfo,
	which int, where Where,
	correctDoubles bool, starts bool) stats.ContingencyTable {
	var actualIn, actualOut int
	var possibleIn, possibleOut int

	for _, pd := range posInfo {
		var in0, in1 bool
		if starts {
			in0 = pd.StartsSite[0]
			in1 = pd.StartsSite[which]
		} else {
			in0 = pd.InSite[0]
			in1 = pd.InSite[which]
		}
		score := findScore(where, correctDoubles, in0, in1)

		if pd.Actual[which] {
			if score > 0 {
				actualIn += score
			} else {
				actualOut++
			}
		}

		if pd.Possible > 0 {
			if score > 0 {
				possibleIn += score * pd.Possible
			} else {
				possibleOut += pd.Possible
			}
		}
	}

	var ret stats.ContingencyTable
	ret.Init(actualIn, actualOut, possibleIn, possibleOut)

	// This gives you a different CT but the same exact OR and p, as it should
	// be:
	// ret.Init(actualIn, possibleIn, actualOut, possibleOut)
	return ret
}

func FindCT(posInfo PosInfo,
	which int, where Where, correctDoubles bool) stats.ContingencyTable {
	return findCT(posInfo, which, where, correctDoubles, false)
}

func FindSiteCT(posInfo PosInfo,
	which int, where Where, correctDoubles bool) stats.ContingencyTable {
	return findCT(posInfo, which, where, correctDoubles, true)
}

// Do it the way they did in the preprint.
func FindCTWrong(posInfo PosInfo,
	which int, where Where, correctDoubles bool) stats.ContingencyTable {
	var silentInside, nsInside int
	var silentOutside, nsOutside int

	for _, pd := range posInfo {
		in0 := pd.InSite[0]
		in1 := pd.InSite[which]
		score := findScore(where, correctDoubles, in0, in1)

		silentMut := pd.Actual[which]

		if silentMut {
			if score > 0 {
				silentInside += score
			} else {
				silentOutside++
			}
		} else {
			if score > 0 {
				nsInside += score
			} else {
				nsOutside++
			}
		}
	}

	var ret stats.ContingencyTable
	ret.Init(silentInside, nsInside, silentOutside, nsOutside)
	return ret
}

type Result struct {
	Genome int
	Sites  [][]byte
	OR     float64
	P      float64
}

/*
Calculate the CT using the default calculation for the BsaI/BsmBI sites
*/
func CalculateCT(g *genomes.Genomes) stats.ContingencyTable {
	calc := CalcCT{FindCT, 1}
	possible := NewPossibleMap(calc.Window,
		mutations.PossibleSilentMuts2(g, 0, calc.Window))
	posInfo := FindPositionInfo(g, possible, RE_SITES)
	return calc.Calc(posInfo, 1, SITE_IN_EITHER, true)
}

/*
Apply num of the mutations. But note that if you used PossibleSilentMuts2 lots
of them won't actually be silent on their own
*/
func Mutate(g* genomes.Genomes,
	possible []mutations.MutatedSequence, num int) {
	n := len(possible)
	done := make(map[int]bool)
	for i := 0; i < num; i++ {
		var j int
		for {
			j = rand.Intn(n)
			if !done[j] {
				break
			}
		}
		mut := possible[j]
		g.Nts[0][mut.Pos] = mut.To[0]
		done[j] = true
	}
}
