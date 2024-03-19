package main

import (
	"fmt"
	"genomics/genomes"
)

func ShowDirections(g *genomes.Genomes,
	trans Transition, concs []Concentration) (int, int) {
	var same, diff int

	for _, conc := range concs {
		a := string(g.Nts[0][conc.Pos:conc.Pos+conc.Length])
		b := string(g.Nts[1][conc.Pos:conc.Pos+conc.Length])

		var direction byte
		if a == trans.ANts && b == trans.BNts {
			direction = 'V'
			same++
		} else if a == trans.BNts && b == trans.ANts {
			direction = '^'
			diff++
		} else {
			continue
		}

		fmt.Printf("%c", direction)
	}
	fmt.Printf(" (%d, %d)\n", same, diff)
	return same, diff
}
