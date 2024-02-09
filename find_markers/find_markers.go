package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"math/rand"
	"os"
	"slices"
)

func findMarkers(g *genomes.Genomes, window int, verbose bool) int {
	var ret int
	nts := g.Nts

positions:
	for i := 0; i < g.Length(); i++ {
		ref := nts[0][i]

		switch ref {
		case 'A':
			fallthrough
		case 'C':
			fallthrough
		case 'G':
			fallthrough
		case 'T':
			break
		default:
			continue positions
		}

		for j := 1; j < g.NumGenomes(); j++ {
			if nts[j][i] == ref {
				continue positions
			}

			if nts[j][i] == '-' {
				continue positions
			}

			// We're looking for all these other g being the same as each
			// other but different from the ref.
			if j > 1 && nts[j][i] != nts[j-1][i] {
				continue positions
			}
		}

		// Now look to see if we are in a conserved region-- is everything the
		// same for some nts either side of this change?
		for j := i - window; j < i+window; j++ {
			if j > 0 && j < g.Length() {
				for k := 2; k < g.NumGenomes(); k++ {
					if nts[k][j] != nts[k-1][j] {
						continue positions
					}
				}
			}
		}

		// We expect the change to be silent.
		if g.HaveOrfs() {
			var env genomes.Environment
			err := env.Init(g, i, 1, 0)
			if err != nil {
				continue
			}
			silent, _ := env.Replace(nts[1][i : i+1])
			if !silent {
				continue
			}
		}

		if verbose {
			fmt.Printf("Position %d: 1st genome has %c, "+
				"the others have %c\n", i, ref, nts[1][i])
		}
		ret++
	}
	return ret
}

func swap(g *genomes.Genomes, i, j int) {
	g.Nts[i], g.Nts[j] = g.Nts[j], g.Nts[i]
	g.Names[i], g.Names[j] = g.Names[j], g.Names[i]
}

func subSample(g *genomes.Genomes, n int) {
	rand.Shuffle(g.NumGenomes(), func(i, j int) {
		swap(g, i, j)
	})

	g.Nts = g.Nts[:n]
	g.Names = g.Names[:n]
}

type Result struct {
	name	string
	count	int
}

type Results map[string]Result

func (r Results) record(name string) {
	record, there := r[name]

	if !there {
		record.name = name
	}

	record.count++
	r[name] = record
}

func (r Results) display() {
	arr := make([]Result, 0, len(r))
	for _, v := range r {
		arr = append(arr, v)
	}
	slices.SortFunc(arr, func(a, b Result) int { return b.count - a.count })

	for _, r := range arr {
		fmt.Printf("%s: %d\n", r.name, r.count)
	}
}

func main() {
	var reorder bool
	var orfs string
	var window int
	var sample int
	var seed int
	var iterations int
	var verbose bool

	flag.BoolVar(&reorder, "r", false, "Try reorderings")
	flag.StringVar(&orfs, "orfs", "", "ORFS file")
	flag.IntVar(&window, "w", 20, "Window for conservedness")
	flag.IntVar(&sample, "s", 0, "Sample this many after the 1st (0 means all)")
	flag.IntVar(&seed, "seed", 1, "Random number seed")
	flag.IntVar(&iterations, "its", 1, "Iterations")
	flag.BoolVar(&verbose, "v", false, "verbose")
	flag.Parse()

	rand.Seed(int64(seed))

	argi := len(os.Args) - flag.NArg()
	g := genomes.LoadGenomes(os.Args[argi], orfs, false)
	g.RemoveGaps()

	results := make(Results)
	orgNts, orgNames := g.Nts, g.Names

	for i := 0; i < iterations; i++ {
		g.Nts, g.Names = orgNts, orgNames
		if sample > 0 {
			subSample(g, sample)
		}

		if verbose {
			for _, name := range g.Names {
				fmt.Println(name)
			}
		}

		n := findMarkers(g, window, verbose)
		if n > 0 {
			results.record(g.Names[0])
		}

		if reorder {
			nts, names := g.Nts, g.Names
			for i := 1; i < g.NumGenomes(); i++ {
				g.Nts, g.Names = nts, names
				swap(g, 0, i)
				if verbose {
					fmt.Printf("%s first\n", g.Names[0])
				}
				n := findMarkers(g, window, verbose)
				if n > 0 {
					results.record(g.Names[0])
				}
			}
		}
	}

	results.display()
}
