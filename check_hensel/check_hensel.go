package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"genomics/utils"
	"log"
	"math/rand"
)

func CalcSimilarities(g *genomes.Genomes, which int) []float64 {
	ret := make([]float64, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret[i] = g.SequenceSimilarity(which, i, false)
	}
	return ret
}

// Return the number of differences. Assume a and b are the same length
func CompareRegion(a []byte, b []byte, start, end int) int {
	aRegion := a[start:end]
	bRegion := b[start:end]

	ret := 0
	for i, val := range aRegion {
		if bRegion[i] != val {
			ret++
		}
	}

	return ret
}

// Describes proximity to a relative
type Proximity struct {
	which       int // Which relative
	differences int // How many differences
	window      int // At what window size
}

// We assume Proximity is sorted by fewest differences first
func HaveUniqueThreeBest(p []Proximity) bool {
	// We just need fresh air between the third and the fourth to satisfy the
	// criterion.
	return p[3].differences > p[2].differences
}

func ShowProximities(g *genomes.Genomes, which int,
	sitePos int, proximities []Proximity) {
	for _, p := range proximities {
		fmt.Println(g.Names[p.which], p.window, p.differences)
	}
}

type Set map[int]bool

// Return the relatives that are closest window either side of the site
func FindClosest(g *genomes.Genomes, which int,
	sitePos int, siteSize int, window int, exclude Set) []Proximity {
	ret := make([]Proximity, 0)

	for i := 0; i < g.NumGenomes(); i++ {
		if i == which {
			continue
		}
		if exclude[i] {
			continue
		}
		end := sitePos
		start := max(0, end-window)
		differences := CompareRegion(g.Nts[which], g.Nts[i], start, end)

		start = min(g.Length(), sitePos+siteSize)
		end = min(g.Length(), start+window)

		differences += CompareRegion(g.Nts[which], g.Nts[i], start, end)
		ret = append(ret, Proximity{i, differences, window})
	}
	utils.SortByKey(ret, true, func(p Proximity) int {
		return p.differences
	})
	return ret
}

type Algo int

const (
	ORIGINAL = iota
	KEEP_BEST
)

// Find the num unambiguous closest either side of the site
func FindNumClosest(g *genomes.Genomes, which int,
	sitePos int, siteSize int, window int, num int, algo Algo) []Proximity {
	num = min(num, g.NumGenomes()-1)
	ret := make([]Proximity, 0, num)
	got := make(Set)

	for {
		need := num - len(ret) // how many more relatives do we need?

		// Find the closest ones, excluding those we already have
		prox := FindClosest(g, which, sitePos, siteSize, window, got)

		// Take any clear leaders
		for _, p := range prox {
			if p.differences < prox[need].differences {
				ret = append(ret, p)
				got[p.which] = true
			}
		}

		if len(ret) == num {
			return ret
		}

		if algo == ORIGINAL {
			/*
				Under the original algorithm, if we didn't find the number we
				wanted (which was three), we threw them all away and then
				widened the window. This meant you threw away some quite good
				matches and ended up with worse ones over a longer window.
			*/
			ret = ret[:0]
			got = make(Set)
		}

		window++
	}
}

// Maps positions to counts of how many of 3 closest relatives matched
type Matches map[int]int

func CompareRelatives(g *genomes.Genomes,
	which int, algo Algo, verbose bool) Matches {
	ret := make(Matches)

	pfn := func(string, ...any) (int, error) {
		return 0, nil
	}
	if verbose {
		pfn = fmt.Printf
	}

	for _, site := range hotspots.RE_SITES {
		for search := genomes.NewLinearSearch(g,
			which, site, 0.0); !search.End(); search.Next() {
			pos, _ := search.Get()
			closest := FindNumClosest(g, which, pos, len(site), 50, 3, algo)
			// Now the question is in how many of the three closest does the
			// actual site also match?
			count := 0
			for _, p := range closest {
				differences := CompareRegion(g.Nts[which], g.Nts[p.which],
					pos, pos+len(site))
				if differences == 0 {
					count++
				}
				pfn("%s at %d in %s (closest over %dnts either "+
					"side) has %d differences in the site\n",
					string(site), pos+1, g.Names[p.which],
					p.window, differences)
			}
			pfn("%d/3 are completely the same\n", count)
			ret[pos] = count
		}
	}
	return ret
}

func AddSite(g *genomes.Genomes, pattern []byte, which int) {
	for {
		pos := rand.Intn(g.Length() - len(pattern))
		for i := pos; i < g.Length()-len(pattern); i++ {
			alreadyThere := true
			for j := 0; j < len(pattern); j++ {
				if g.Nts[which][i+j] != pattern[j] {
					alreadyThere = false
				}
			}
			if alreadyThere {
				continue
			}
			silent, numMuts, _ := genomes.IsSilentWithReplacement(g, i,
				which, which, pattern)
			if silent && numMuts == 1 {
				fmt.Printf("New %s site at %d\n", string(pattern), i)
				copy(g.Nts[which][i:i+len(pattern)], pattern)
				// g.Save("Modified", "modified.fasta", 0)
				return
			}
		}
	}
}

func AddsSite(g *genomes.Genomes,
	which int, m *mutations.Mutation) ([]byte, int) {
	start := max(0, m.Pos-5)
	end := min(start+6, g.Length())
	for i := start; i < end; i++ {
		for _, site := range hotspots.RE_SITES {
			matched := true
			for j := 0; j < len(site); j++ {
				var nt byte
				if i+j == m.Pos {
					nt = m.To
				} else {
					nt = g.Nts[which][i+j]
				}
				if nt != site[j] {
					matched = false
				}
			}
			if matched {
				return site, i
			}
		}
	}
	return nil, 0
}

func Estimate(g *genomes.Genomes, algo Algo) {
	fmt.Println("All possible silent point muts that add a site:")
	numPossible, numAdd, numMin1, numMin2, numMin3 := 0, 0, 0, 0, 0
	for which := 0; which < g.NumGenomes(); which++ {
		possible := mutations.PossibleSilentMuts(g, which)
		numPossible += len(possible)

		for _, mut := range possible {
			if site, pos := AddsSite(g, which, &mut); site != nil {
				numAdd++
				// Apply the mutation
				g.Nts[which][mut.Pos] = mut.To
				matches := CompareRelatives(g, which, algo, false)[pos]
				if matches != 0 {
					numMin1++
					fmt.Printf("%c%d%c adds site %s into %s at %d"+
						": %d matches in relatives\n",
						mut.From, mut.Pos+1, mut.To, string(site),
						g.Names[which], pos+1, matches)
				}
				if matches >= 2 {
					numMin2++
				}
				if matches >= 3 {
					numMin3++
				}
				// Unapply point mutation again
				g.Nts[which][mut.Pos] = mut.From
			}
		}
	}
	fmt.Printf("%d possible silent mutations, "+
		"%d of which add sites of which %d match relatives "+
		"(min 1) or %d (min 2) or %d (min 3)\n",
		numPossible, numAdd, numMin1, numMin2, numMin3)
}

func main() {
	var (
		tamper   bool
		estimate bool
		which    int
		algoS    string
		algo     Algo
		fasta    string
		orfs     string
	)

	flag.BoolVar(&tamper, "tamper", false, "Whether to adjust sites")
	flag.BoolVar(&estimate, "estimate", false, "Estimate out of possible muts")
	flag.IntVar(&which, "which", 0, "Which genome to examine")
	flag.StringVar(&algoS, "algo", "keep_best", "Algorithm")
	flag.StringVar(&fasta, "fasta",
		"../fasta/Hassanin.fasta", "Alignment to use")
	flag.StringVar(&orfs, "Orfs", "../fasta/WH1.orfs", "Orfs")
	flag.Parse()

	switch algoS {
	case "original":
		algo = ORIGINAL
	case "keep_best":
		algo = KEEP_BEST
	default:
		log.Fatal("Don't know that algo")
	}

	g := genomes.LoadGenomes(fasta, orfs, false)

	if estimate {
		Estimate(g, algo)
	}

	if tamper {
		for _, site := range hotspots.RE_SITES {
			AddSite(g, site, which)
		}
	}

	fmt.Printf("Comparing sites in %s\n", g.Names[which])
	CompareRelatives(g, which, algo, true)
}
