package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/simulation"
	"genomics/stats"
	"genomics/utils"
	"log"
	"math/rand"
	"os"
)

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

/*
Calculate a Contingency Table in a variety of correct or incorrect ways given
the actual mutations, the possible mutations, where the sites are, and the
length of the genome. If correctDoubles then count where sites appear in both
genomes correctly (only applies when where == SITE_IN_EITHER)
*/
type CalcCT func(posInfo PosInfo,
	which int, where Where, correctDoubles bool) stats.ContingencyTable

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
			score++
		}
	}
	return score
}

func FindCT(posInfo PosInfo,
	which int, where Where, correctDoubles bool) stats.ContingencyTable {
	var actualIn, actualOut int
	var possibleIn, possibleOut int

	for _, pd := range posInfo {
		in0 := pd.InSite[0]
		in1 := pd.InSite[which]
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
				possibleIn += score
			} else {
				possibleOut++
			}
		}
	}

	var ret stats.ContingencyTable
	ret.Init(actualIn, actualOut, possibleIn, possibleOut)
	return ret
}

func FindCTAlt(posInfo PosInfo,
	which int, correctDoubles bool) stats.ContingencyTable {
	var actualIn, actualOut int
	var possibleIn, possibleOut int

	for _, pd := range posInfo {
		in0 := pd.InSite[0]
		in1 := pd.InSite[which]
		actual := pd.Actual[which]

		// Another way of writing the same code (with correctDoubles) which is
		// equivalent but perhaps clearer?
		if actual {
			if in0 {
				actualIn++
			}
			if in1 {
				actualIn++
			}
			if !(in0 || in1) {
				actualOut++
			}
		}

		if pd.Possible > 0 {
			if in0 {
				possibleIn++
			}
			if in1 {
				possibleIn++
			}
			if !(in0 || in1) {
				possibleOut++
			}
		}
	}

	var ret stats.ContingencyTable
	ret.Init(actualIn, actualOut, possibleIn, possibleOut)
	return ret
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

/*
Find a CT based on the number of sites that have been altered, regardless of
how many mutations per site
*/
func FindSiteCT(posInfo PosInfo,
	which int, where Where, correctDoubles bool) stats.ContingencyTable {
	var actualIn, actualOut int
	var possibleIn, possibleOut int

	for i := 0; i < len(posInfo)-6; i++ {
		pd := posInfo[i]

		in0 := pd.StartsSite[0]
		in1 := pd.StartsSite[which]
		score := findScore(where, correctDoubles, in0, in1)

		// Are there >1 actual muts in the whole site?
		var got bool
		for j := 0; j < 6; j++ {
			if posInfo[i+j].Actual[which] {
				got = true
				break
			}
		}
		if got {
			if score > 0 {
				actualIn += score
			} else {
				actualOut++
			}
		}

		// Same thing for possible
		got = false
		for j := 0; j < 6; j++ {
			if posInfo[i+j].Possible > 0 {
				got = true
				break
			}
		}
		if got {
			if score > 0 {
				possibleIn += score
			} else {
				possibleOut++
			}
		}
	}

	var ret stats.ContingencyTable
	ret.Init(actualIn, actualOut, possibleIn, possibleOut)
	return ret
}

type Result struct {
	genome int
	OR     float64
	p      float64
}

/*
If useSites use random sites, and find where they are. Otherwise just use
random positions
*/
func MonteCarlo(g *genomes.Genomes, possible PossibleMap,
	its int, calc CalcCT, where Where, correctDoubles bool) []Result {
	ret := make([]Result, 0)
	for i := 0; i < its; i++ {
		sites := make([][]byte, 4)
		for j := 0; j < 2; j++ {
			sites[j] = utils.RandomNts(6)
			sites[j+2] = utils.ReverseComplement(sites[j])
		}
		posInfo := FindPositionInfo(g, possible, sites)

		for j := 1; j < g.NumGenomes(); j++ {
			ct := calc(posInfo, j, where, correctDoubles)
			OR, p := ct.FisherExact(stats.GREATER)
			ret = append(ret, Result{j, OR, p})
		}

		if i%100 == 0 {
			fmt.Printf("%d/%d\n", i, its)
		}
	}
	return ret
}

func OutputResults(results []Result, which int) {
	f, err := os.Create("ORs")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	orW := bufio.NewWriter(f)

	f2, err := os.Create("ps")
	if err != nil {
		log.Fatal(err)
	}
	defer f2.Close()
	pW := bufio.NewWriter(f2)

	for _, res := range results {
		if res.genome == which {
			fmt.Fprintln(orW, res.OR)
			fmt.Fprintln(pW, res.p)
		}
	}

	orW.Flush()
	pW.Flush()

	fmt.Println("Wrote ORs and ps")
}

func MakeTestGenomes(g *genomes.Genomes) {
	nd := mutations.NewNucDistro(g)

	var orfs genomes.Orfs
	ret := genomes.NewGenomes(orfs, 2)

	ret.Nts[0] = make([]byte, g.Length())
	ret.Nts[1] = make([]byte, g.Length())

	for i := 0; i < g.Length(); i++ {
		ret.Nts[0][i] = nd.Random()
		if rand.Float64() < 0.96 {
			ret.Nts[1][i] = ret.Nts[0][i]
		} else {
			ret.Nts[1][i] = nd.Random()
		}
	}
	ret.Names[0] = "Test Random Genome"
	ret.Names[1] = "Test Random Relative"

	ret.SaveMulti("random.fasta")
}

func main() {
	sites := [][]byte{
		[]byte("GGTCTC"),
		[]byte("GAGACC"),
		[]byte("CGTCTC"),
		[]byte("GAGACG"),
	}

	var fasta string
	var orfs string
	var doMC bool
	var correctDoubles bool
	var algorithm string
	var whichMC int
	var redistribute bool
	var whereS string
	var show bool

	flag.StringVar(&fasta, "fasta",
		"../fasta/CloseRelatives.fasta", "relatives")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "orfs")
	flag.BoolVar(&doMC, "montecarlo", false, "Do the MonteCarlo")
	flag.BoolVar(&correctDoubles, "doubles", true, "Count doubles correctly")
	flag.StringVar(&algorithm, "algo", "default", "Algorithm for finding CT")
	flag.IntVar(&whichMC, "which-mc", 2, "Which genome to use for MonteCarlo")
	flag.BoolVar(&redistribute, "redistrib", false, "Redistribute mutations")
	flag.StringVar(&whereS, "where", "either", "Where to look for sites")
	flag.BoolVar(&show, "show", false, "Show the info")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	if redistribute {
		fmt.Println("Redistributing the mutations")
		nd := mutations.NewNucDistro(g)
		g, _ = simulation.MakeSimulatedMutant(g, 0, whichMC, nd)
		whichMC = 1
	}

	var where Where
	switch whereS {
	case "either":
		where = SITE_IN_EITHER
	case "a":
		where = SITE_IN_A
	case "b":
		where = SITE_IN_B
	default:
		log.Fatal("Invalid where specification")
	}

	possible := NewPossibleMap(mutations.PossibleSilentMuts(g, 0))
	posInfo := FindPositionInfo(g, possible, sites)

	if show {
		posInfo.SaveTSV()
		posInfo.ShowSites()
		return
	}

	var calcFn CalcCT
	switch algorithm {
	case "default":
		// calcFn = FindCTAlt
		calcFn = FindCT
	case "per-site":
		calcFn = FindSiteCT
	case "wrong":
		calcFn = FindCTWrong
	default:
		log.Fatal("Bad algorithm")
	}

	for i := 1; i < g.NumGenomes(); i++ {
		ct := calcFn(posInfo, i, where, correctDoubles)
		OR, p := ct.FisherExact(stats.GREATER)
		fmt.Printf("%s: %f %f\n", g.Names[i], OR, p)
		fmt.Println(ct.String())
	}

	if doMC {
		g = g.Filter(0, whichMC)
		mc := MonteCarlo(g, possible, 5000, calcFn, where, correctDoubles)
		OutputResults(mc, 1)
	}
}
