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
	name  string
	which int // which genome from the original set this is
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

func (p *Pattern) NumMuts(other *Pattern) int {
	var ret int
	for i := 0; i < len(p.nts); i++ {
		if p.nts[i] != other.nts[i] {
			ret++
		}
	}
	return ret
}

func ByLocation(r Results) map[int][]Pattern {
	ret := make(map[int][]Pattern)
	for _, s := range r {
		for k, _ := range s {
			ret[k.pos] = append(ret[k.pos], k)
		}
	}
	return ret
}

func mapPattern(g *genomes.Genomes, start, length int, results Results) {
	// Maps pattern to how many genomes have that pattern
	m := make(map[string]int)

	// Maps pattern to the indices of the genomes with that pattern at this
	// location.
	indices := make(map[string][]int)

	for i := 0; i < g.NumGenomes(); i++ {
		pat := string(g.Nts[i][start : start+length])
		m[pat]++
		indices[pat] = append(indices[pat], i)
	}

	for k, v := range m {
		if !utils.IsRegularPattern([]byte(k)) {
			continue
		}
		for _, index := range indices[k] {
			name := g.Names[index]
			ps, there := results[name]
			if !there {
				ps = make(PatternSet)
			}
			ps[Pattern{name, index, k, start, v}] = true
			results[name] = ps
		}
	}
}

/*
Look for anywhere any of the patterns appears in an alignment, or a
sequence that could be silently mutated into one, and see if any of those
are unique across the set
*/
func FindPatterns(g *genomes.Genomes, patterns []string) Results {
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

			mapPattern(g, j, 6, ret)
		}
	}
	return ret
}

func mostDifferent(patterns []Pattern) (*Pattern, *Pattern, int) {
	mostMuts := 0
	var bestP, bestQ *Pattern

	for i, p := range patterns {
		for j, q := range patterns {
			nm := p.NumMuts(&q)
			if nm > mostMuts {
				mostMuts = nm
				bestP, bestQ = &patterns[i], &patterns[j]
			}
		}
	}
	return bestP, bestQ, mostMuts
}

func main() {
	/*
		g := genomes.LoadGenomes("../fasta/more_relatives.fasta",
			"../fasta/WH1.orfs", false)
	*/
	/*
	g := genomes.LoadGenomes("../fasta/relatives2.fasta",
		"../fasta/WH1.orfs", false)
	*/
	g := genomes.LoadGenomes("../fasta/more_relatives2.fasta",
		"../fasta/WH1.orfs", false)

	interesting := []string{
		"CGTCTC",
		"GAGACC",
		"GGTCTC",
		"GAGACG",
	}

	results := FindPatterns(g, interesting)
	//results.Print()

	interestingSet := ToSet(interesting)
	byLocation := ByLocation(results)
	for k, v := range byLocation {
		p, q, md := mostDifferent(v)

		for _, pat := range v {
			fmt.Printf("%d: %s has %s\n", k, pat.name, pat.nts)
		}

		if md >= 4 && (interestingSet[p.nts] || interestingSet[q.nts]) {
			ss := g.SequenceSimilarity(p.which, q.which) * 100
			fmt.Printf("%d %s vs %s: %s/%s %d muts %.2f%% ss\n", k,
				p.name, q.name, p.nts, q.nts, md, ss)
		}
	}
}
