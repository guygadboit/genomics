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

func checkUnique(g *genomes.Genomes, start, length int) {
	// Maps pattern to how many genomes have that pattern
	m := make(map[string]int)

	// Maps pattern to the last seen name with that pattern. We don't bother
	// with a list of names because we're only interested in the unique ones.
	names := make(map[string]string)

	for i := 0; i < g.NumGenomes(); i++ {
		pat := string(g.Nts[i][start:start+length])
		m[pat]++
		names[pat] = g.Names[i]
	}

	for k, v := range m {
		if !utils.IsRegularPattern([]byte(k)) {
			continue
		}
		fmt.Printf("%s at %d in <%s> found %d times\n", k, start, names[k], v)
	}
}

/*
	Look for anywhere any of the patterns appears in an alignment, or a
	sequence that could be silently mutated into one, and see if any of those
	are unique across the set
*/
func FindUnique(g *genomes.Genomes, patterns []string) {
	patSet := ToSet(patterns)

	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < g.Length(); j++ {
			var env genomes.Environment
			err := env.Init(g, j, 6, i)
			if err != nil {
				continue
			}

			// Do we actually have the pattern here?
			sp := string(g.Nts[i][j:j+6])
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

			checkUnique(g, j, 6)
		}
	}
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

	FindUnique(g, interesting)
}
