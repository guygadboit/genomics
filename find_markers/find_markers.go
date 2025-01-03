package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"math/rand"
	"os"
	"slices"
)

func analyzeIsland(g *genomes.Genomes, pos int) {
	counter := make(map[byte]int)
	for i := 0; i < g.NumGenomes(); i++ {
		nt := g.Nts[i][pos]
		count, _ := counter[nt]
		counter[nt] = count + 1
	}
	for k, v := range counter {
		fmt.Printf("%s: %d\n", string(k), v)
	}
}

type findFunc func(g *genomes.Genomes, window int, verbose bool) int

func isRegularNt(nt byte) bool {
	switch nt {
	case 'A':
		fallthrough
	case 'C':
		fallthrough
	case 'G':
		fallthrough
	case 'T':
		return true
	default:
		return false
	}
}

// Look for SNPs in conserved areas
func findIslands(g *genomes.Genomes, window int, verbose bool) int {
	nts := g.Nts
	var ret int
	muts := make([]bool, g.Length())

	for i := 0; i < g.Length(); i++ {
		for j := 1; j < g.NumGenomes(); j++ {
			a, b := nts[j][0], nts[j][i]
			if isRegularNt(a) && isRegularNt(b) {
				if a != b {
					muts[i] = true
					break
				}
			}
		}
	}

	var conserved int
	for i, mut := range muts {
		if mut && conserved >= window {
			// Now check for window the other side
			good := true
			for j := i + 1; j < i+window; j++ {
				if j >= g.Length() {
					good = false
					break
				}
				if muts[j] {
					good = false
					break
				}
			}
			if good {
				ret++
				fmt.Printf("Island at %d\n", i)
				analyzeIsland(g, i)
			}
		}

		if mut {
			conserved = 0
		} else {
			conserved++
		}
	}
	return ret
}

/*
This is looking for something that is different in the first genome but the
same in all the others
*/
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

			// We're looking for all these others being the same as each other
			// but different from the ref.
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
			fmt.Printf("Position %d: %s has %c, "+
				"the others have %c\n", i, g.Names[0], ref, nts[1][i])
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
	name  string
	count int
}

type Results map[string]Result

func (r Results) record(name string, nts int) {
	record, there := r[name]

	if !there {
		record.name = name
	}

	record.count += nts
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

/*
We keep subsampling and shuffling genomes around. So maintain a stack so we can
keep reverting those operations back.
*/
type Stack struct {
	data []genomes.Genomes
}

func (s *Stack) push(g *genomes.Genomes) {
	s.data = append(s.data, *g) // note we're taking a copy
}

func (s *Stack) pop() *genomes.Genomes {
	n := len(s.data) - 1
	ret := s.data[n]
	s.data = s.data[:n]
	return &ret
}

func averageSimilarity(g *genomes.Genomes) float64 {
	var total float64
	for i := 1; i < g.NumGenomes(); i++ {
		total += g.SequenceSimilarity(0, i)
	}
	return total / float64(g.NumGenomes()-1)
}

func main() {
	var reorder bool
	var orfs string
	var window int
	var sample int
	var seed int
	var iterations int
	var verbose bool
	var islands bool

	flag.BoolVar(&reorder, "r", false, "Try reorderings")
	flag.StringVar(&orfs, "orfs", "", "ORFS file")
	flag.IntVar(&window, "w", 0, "Window for conservedness (0 means auto)")
	flag.IntVar(&sample, "s", 0, "Sample this many after the 1st (0 means all)")
	flag.IntVar(&seed, "seed", 1, "Random number seed")
	flag.IntVar(&iterations, "its", 1, "Iterations")
	flag.BoolVar(&verbose, "v", false, "verbose")
	flag.BoolVar(&islands, "islands", false, "Look for islands")
	flag.Parse()

	rand.Seed(int64(seed))

	var stack Stack
	argi := len(os.Args) - flag.NArg()
	g := genomes.LoadGenomes(os.Args[argi], orfs, false)
	g.RemoveGaps()

	results := make(Results)

	var fn findFunc
	if islands {
		fn = findIslands
	} else {
		fn = findMarkers
	}

	for i := 0; i < iterations; i++ {
		stack.push(g)
		if sample > 0 {
			subSample(g, sample)
		}

		// Set the window based on the similarity-- we're looking for isolated
		// muts in big seas of conservedness
		if window == 0 {
			similarity := averageSimilarity(g)
			window = int(50 / similarity)
		}

		/*
			if verbose {
				for _, name := range g.Names {
					fmt.Println(name)
				}
			}
		*/

		n := fn(g, window, verbose)
		if n > 0 {
			results.record(g.Names[0], n)
		}

		if reorder {
			for i := 1; i < g.NumGenomes(); i++ {
				stack.push(g)
				swap(g, 0, i)
				n := fn(g, window, verbose)
				if n > 0 {
					results.record(g.Names[0], n)
				}
				g = stack.pop()
			}
		}
		g = stack.pop()
	}

	results.display()
}
