package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"os"
)

func findMarkers(g *genomes.Genomes) {
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
		for j := i - 20; j < i+20; j++ {
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

		fmt.Printf("Position %d: 1st genome has %c, "+
			"the others have %c\n", i, ref, nts[1][i])
	}
}

func swap(g *genomes.Genomes, i, j int) {
	g.Nts[i], g.Nts[j] = g.Nts[j], g.Nts[i]
	g.Names[i], g.Names[j] = g.Names[j], g.Names[i]
}

func main() {
	var reorder bool
	var orfs string

	flag.BoolVar(&reorder, "r", false, "Try reorderings")
	flag.StringVar(&orfs, "orfs", "", "ORFS file")
	flag.Parse()

	argi := len(os.Args) - flag.NArg()
	g := genomes.LoadGenomes(os.Args[argi], orfs, false)
	g.RemoveGaps()

	for _, name := range g.Names {
		fmt.Println(name)
	}

	findMarkers(g)

	if reorder {
		orgNts := g.Nts
		for i := 1; i < g.NumGenomes(); i++ {
			g.Nts = orgNts
			swap(g, 0, i)
			fmt.Printf("%s first\n", g.Names[0])
			findMarkers(g)
		}
	}
}
