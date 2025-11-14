package main

import (
	"flag"
	"fmt"
	. "genomics/outgroup"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/mutations"
	"log"
	"math/rand"
)

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
				differences := CompareRegion(g.Nts[which], g.Nts[p.Which],
					pos, pos+len(site))
				if differences == 0 {
					count++
				}
				pfn("%s at %d in %s (closest over %dnts either "+
					"side) has %d differences in the site\n",
					string(site), pos+1, g.Names[p.Which],
					p.Window, differences)
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
