package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

type StringSet map[string]bool

func ToSet(st []string) StringSet {
	ret := make(StringSet)
	for _, s := range st {
		ret[s] = true
	}
	return ret
}

type Pattern struct {
	nts   string
	pos   int
	count int // how many relatives it is in
}

type PatternSet map[Pattern]bool

// Maps each genome to its patterns
type Results map[string]PatternSet

func (r Results) Print() {
	for k, v := range r {
		fmt.Printf("%s\n", k)
		/*
		slices.SortFunc(v, func(a, b Pattern) int {
			return a.count - b.count
		})
		*/
		for k, _ := range v {
			fmt.Printf("%s %d %d\n", k.nts, k.pos, k.count)
		}
	}
}

func checkUnique(g *genomes.Genomes, start, length int, results Results) {

	// Maps pattern to how many genomes have that pattern
	m := make(map[string]int)

	// Maps pattern to the names of the genomes with that pattern at this
	// location.
	names := make(map[string][]string)

	for i := 0; i < g.NumGenomes(); i++ {
		pat := string(g.Nts[i][start : start+length])
		m[pat]++
		names[pat] = append(names[pat], g.Names[i])
	}

	for k, v := range m {
		if !utils.IsRegularPattern([]byte(k)) {
			continue
		}
		for _, name := range names[k] {
			ps, there := results[name]
			if !there {
				ps = make(PatternSet)
			}
			ps[Pattern{k, start, v}] = true
			results[name] = ps
		}
	}
}

/*
Look for anywhere any of the patterns appears in an alignment, or a
sequence that could be silently mutated into one, and see if any of those
are unique across the set
*/
func FindUnique(g *genomes.Genomes, patterns []string) Results {
	ret := make(Results)

	patSet := ToSet(patterns)

	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < g.Length(); j++ {
			var env genomes.Environment
			err := env.Init(g, j, 6, i)
			if err != nil {
				continue
			}

			// Do we actually have the pattern here?
			sp := string(g.Nts[i][j : j+6])
			interested := patSet[sp]

			// Or is one of the silent alternatives the pattern?
			if !interested {
				for _, alt := range env.FindAlternatives(6) {
					if patSet[string(alt.Nts)] {
						interested = true
						break
					}
				}
			}

			if !interested {
				continue
			}

			checkUnique(g, j, 6, ret)
		}
	}
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/more_relatives.fasta",
		"../fasta/WH1.orfs", false)

	interesting := []string{
		"CGTCTC",
		"GAGACC",
		"GGTCTC",
		"GAGACG",
	}

	results := FindUnique(g, interesting)
	results.Print()
}
