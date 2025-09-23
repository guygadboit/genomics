package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"strings"
)

var TrsMarkers string
var RbdMarkers string

func LoadTrsMarkers() {
	utils.Lines("./trs_locations", func(line string, err error) bool {
		location := utils.Atoi(line) - 1
		// Those are all the places where ACGAAC occurs, which is the TRS
		// marker for SARS-CoV-2 (and also SARS1, and I think all of these).
		TrsMarkers = fmt.Sprintf(`%s
set arrow from %d,graph 0 to %d,graph 1 nohead filled lc "red"`,
			TrsMarkers, location, location)
		return true
	})
}

func LoadRbdMarkers() {
	// Note these are only applicable to a SARS2 alignment!
	/*
			RbdMarkers = fmt.Sprintf(`
		set arrow from 22874, graph 0 to 22874, graph 1 nohead filled lc "red"
		set arrow from 23081, graph 0 to 23081, graph 1 nohead filled lc "red"
		set arrow from 22558, graph 0 to 22558, graph 1 nohead filled lc "blue"
		set arrow from 23600, graph 0 to 23600, graph 1 nohead filled lc "blue"
		`)
	*/

	/*
		SARS1 RBM
		set arrow from 22763, graph 0 to 22763, graph 1 nohead filled lc "blue"
		set arrow from 22968, graph 0 to 22968, graph 1 nohead filled lc "blue"
	*/

	// SARS1 RBD
	RbdMarkers = fmt.Sprintf(`
set arrow from 22448, graph 0 to 22448, graph 1 nohead filled lc "red"
set arrow from 23487, graph 0 to 23487, graph 1 nohead filled lc "red"
`)
}

func ORFMarkers(g *genomes.Genomes, which ...string) string {
	include := utils.ToSet(which)
	ret := make([]string, 0, 2*len(which))
	colours := []string{"red", "blue", "green"}

	for i, orf := range g.Orfs {
		if include[orf.Name] {
			colour := colours[i%len(colours)]
			ret = append(ret, fmt.Sprintf("set arrow from %d, graph 0 to %d,"+
				" graph 1 nohead filled lc \"%s\"",
				orf.Start, orf.Start, colour))
			ret = append(ret, fmt.Sprintf("set arrow from %d, graph 0 to %d,"+
				" graph 1 nohead filled lc \"%s\"",
				orf.End, orf.End, colour))
		}
	}
	return strings.Join(ret, "\n")
}

func JumpMarker(start, end int) string {
	return fmt.Sprintf(`
set arrow from %d, graph 0 to %d, graph 1 nohead filled lc "magenta"
set arrow from %d, graph 0 to %d, graph 1 nohead filled lc "magenta"`,
		start, start, end, end)
}

func init() {
	LoadTrsMarkers()
	LoadRbdMarkers()
}
