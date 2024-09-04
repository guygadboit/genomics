package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	//"genomics/simulation"
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

type CalcCT struct {
	/*
		Calculate a Contingency Table in a variety of correct or incorrect ways
		given the actual mutations, the possible mutations, where the sites are,
		and the length of the genome. If correctDoubles then count where sites
		appear in both genomes correctly (only applies when where ==
		SITE_IN_EITHER)
	*/
	calc func(posInfo PosInfo, which int,
		where Where, correctDoubles bool) stats.ContingencyTable

	// How long a window to use when looking at actually and possibly silently
	// mutated sequences.
	window int
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
				// possibleIn += score
				possibleIn += score * pd.Possible
			} else {
				// possibleOut++
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
	genome int
	sites  [][]byte
	OR     float64
	p      float64
}

/*
If useSites use random sites, and find where they are. Otherwise just use
random positions
*/
func MonteCarlo(g *genomes.Genomes, possible *PossibleMap,
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
			ct := calc.calc(posInfo, j, where, correctDoubles)
			OR, p := ct.FisherExact(stats.GREATER)
			ret = append(ret, Result{j, sites, OR, p})
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
			if res.OR != 0 {
				fmt.Fprintln(orW, res.OR)
			}
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

// Redistribute the silent mutations randomly according to position. Do this in
// all the genomes.
func Redistribute(g *genomes.Genomes, possible *PossibleMap) *genomes.Genomes {
	ret := g.Clone()

	if possible.window != 1 {
		log.Fatal("Can't do this with non-1 based possible maps")
	}

	for i := 1; i < ret.NumGenomes(); i++ {
		g2 := ret.Filter(0, i)
		numSilent, _ := mutations.CountMutations(g2)
		fmt.Printf("There are %d silent muts\n", numSilent)

		positions := make([]int, 0, len(possible.mutations))
		for k, _ := range possible.mutations {
			positions = append(positions, k)
		}

		// They're in a fairly random order anyway (because map keys) but
		// shuffle them again to be sure.
		rand.Shuffle(len(positions), func(i, j int) {
			positions[i], positions[j] = positions[j], positions[i]
		})

		// Set the i'th genome to a copy of the 0th, ready to apply the random
		// mutations.
		ret.Nts[i] = ret.Nts[0]
		ret.DeepCopy(i)

		mutsToApply := numSilent
		for j := 0; j < len(positions); j++ {
			muts := possible.mutations[positions[j]]
			k := rand.Intn(len(muts))
			mut := muts[k]
			ret.Nts[i][mut.Pos] = mut.To[0]
			mutsToApply--
			if mutsToApply == 0 {
				break
			}
		}

		if mutsToApply != 0 {
			// This is pretty unlikely
			fmt.Printf("Wasn't able to apply all of them. %d left\n",
				mutsToApply)
		}
	}

	return ret
}

// Highlight the sites
func makeHighlights(pi PosInfo, which int) []genomes.Highlight {
	ret := make([]genomes.Highlight, 0)

	for pos, pd := range pi {
		in0 := pd.InSite[0]
		in1 := pd.InSite[which]

		if in0 && in1 {
			ret = append(ret, genomes.Highlight{pos, pos + 1, 'x'})
		} else if in0 {
			ret = append(ret, genomes.Highlight{pos, pos + 1, 'w'})
		} else if in1 {
			ret = append(ret, genomes.Highlight{pos, pos + 1, 'v'})
		}
	}

	return ret
}

func TestAll(g *genomes.Genomes, sites [][]byte,
	calc CalcCT, where Where, correctDoubles bool, redistribute bool) []Result {
	count := 0
	ret := make([]Result, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		g2 := g.Swap(0, i)
		a := i
		possible := NewPossibleMap(calc.window,
			mutations.PossibleSilentMuts2(g2, 0, calc.window))
		posInfo := FindPositionInfo(g2, possible, sites)
		for j := 1; j < g2.NumGenomes(); j++ {
			b := j
			if j == i {
				b = 0
			}
			ct := calc.calc(posInfo, j, where, correctDoubles)
			OR := ct.CalcOR()
			fmt.Printf("%d,%d: %f", a, b, OR)
			if ct.OR > 3 {
				_, p := ct.FisherExact(stats.GREATER)
				fmt.Printf(" p=%f\n", p)
			} else {
				fmt.Printf("\n")
			}
			count++

			// We aren't currently working out the p for all of them, and we're
			// just putting 0 in as the genome. If you ever need to do anything
			// except output all the ORs change this.
			ret = append(ret, Result{0, sites, ct.OR, ct.P})
		}
	}
	fmt.Printf("Tested %d pairs (1/%d=%f)\n", count, count, 1.0/float64(count))
	return ret
}

func main() {
	sites := [][]byte{
		[]byte("GGTCTC"),
		[]byte("GAGACC"),
		[]byte("CGTCTC"),
		[]byte("GAGACG"),
		/*
			[]byte("CGTCTT"),
			[]byte("AAGACG"),
			[]byte("TTGTTA"),
			[]byte("TAACAA"),
		*/
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
	var its int
	var testAll bool
	var save bool

	flag.StringVar(&fasta, "fasta",
		"../fasta/CloseRelatives.fasta", "relatives")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "orfs")
	flag.BoolVar(&doMC, "montecarlo", false, "Do the MonteCarlo")
	flag.IntVar(&its, "its", 5000, "Montecarlo iterations")
	flag.BoolVar(&correctDoubles, "doubles", true, "Count doubles correctly")
	flag.StringVar(&algorithm, "algo", "default", "Algorithm for finding CT")
	flag.IntVar(&whichMC, "which-mc", 2, "Which genome to use for MonteCarlo")
	flag.BoolVar(&redistribute, "redistrib", false, "Redistribute mutations")
	flag.StringVar(&whereS, "where", "either", "Where to look for sites")
	flag.BoolVar(&show, "show", false, "Show the info")
	flag.BoolVar(&testAll, "testall", false, "Test all pairs")
	flag.BoolVar(&save, "save", false, "Save clu files")

	flag.Parse()

	var calc CalcCT
	switch algorithm {
	case "default":
		// calcFn = FindCTAlt
		calc = CalcCT{FindCT, 1}
	case "per-site":
		calc = CalcCT{FindSiteCT, 6}
	case "wrong":
		calc = CalcCT{FindCTWrong, 1}
	default:
		log.Fatal("Bad algorithm")
	}

	g := genomes.LoadGenomes(fasta, orfs, false)
	possible := NewPossibleMap(calc.window,
		mutations.PossibleSilentMuts2(g, 0, calc.window))

	if redistribute {
		fmt.Println("Redistributing the mutations")
		muts := mutations.ToSequences(mutations.PossibleSilentMuts(g, 0))
		pm := NewPossibleMap(1, muts)
		g = Redistribute(g, pm)
		g.SaveMulti("redistributed.fasta")
		fmt.Printf("Wrote redistributed.fasta\n")
	}

	posInfo := FindPositionInfo(g, possible, sites)

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

	if show {
		posInfo.SaveTSV()
		posInfo.ShowSites(g)
	}

	if doMC {
		g2 := g.Filter(0, whichMC)
		mc := MonteCarlo(g2, possible, its, calc, where, correctDoubles)
		OutputResults(mc, 1)

		var greater int
		ct := calc.calc(posInfo, whichMC, where, correctDoubles)
		refOR := ct.CalcOR()
		for _, result := range mc {
			if result.OR >= refOR {
				greater++
				fmt.Printf("%s/%s %f\n",
					string(result.sites[0]),
					string(result.sites[1]), result.OR)
			}
		}
		fmt.Printf("%f are greater than reference OR of %f\n",
			float64(greater)/float64(its), refOR)
	}

	for i := 1; i < g.NumGenomes(); i++ {
		ct := calc.calc(posInfo, i, where, correctDoubles)
		OR, p := ct.FisherExact(stats.GREATER)
		fmt.Printf("%s: %f %.4g\n", g.Names[i], OR, p)
		fmt.Println(ct.HumanString())

		if save {
			highlights := makeHighlights(posInfo, i)
			fname := fmt.Sprintf("%d.clu", i)
			g.SaveWithTranslation(fname, highlights, 0, i)
		}
	}

	if testAll {
		results := TestAll(g, sites, calc, where, correctDoubles, false)
		OutputResults(results, 0)
	}
}
