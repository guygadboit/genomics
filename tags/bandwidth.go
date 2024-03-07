package main

import (
	"fmt"
	"genomics/genomes"
)

// Return all possible combinations of window nts
func AllPatterns(window int) [][]byte {
	nts := []byte("GATC")
	ret := make([][]byte, 0)
	current := make([]byte, window)

	toPat := func(inp []byte) []byte {
		ret := make([]byte, len(inp))
		for i, v := range inp {
			ret[i] = nts[v]
		}
		return ret
	}

all:
	for {
		ret = append(ret, toPat(current))

		for i := 0; i < len(current); i++ {
			current[i]++
			if current[i] < 4 {
				break
			}
			if i == window-1 {
				break all
			}
			for j := 0; j <= i; j++ {
				current[j] = 0
			}
		}
	}
	return ret
}

type NtRange struct {
	window int
	nts    [][]byte
}

func (nr *NtRange) Init(window int) {
	nr.window = window
	nr.nts = AllPatterns(window)
}

// How many ways are there of changing window nts at pos silently? Return the
// maximum number of muts.
func (nr *NtRange) Alternatives(g *genomes.Genomes,
	which int, pos int) int {
	var env genomes.Environment

	err := env.Init(g, pos, nr.window, which)
	if err != nil {
		return 0
	}

	max := 0

	for _, alt := range nr.nts {
		silent, numMuts := env.Replace(alt)
		if silent {
			if numMuts > max {
				max = numMuts
			}
		}
	}

	return max
}

// The maximum numbers of muts available at each position
func Bandwidth(g *genomes.Genomes, which int,
	window int, verbose bool) []int {
	ret := make([]int, g.Length())

	var nr NtRange
	nr.Init(window)

	for i := 0; i < g.Length(); i++ {
		ret[i] = nr.Alternatives(g, which, i)
		if verbose && i%1000 == 0 {
			fmt.Printf("%d / %d\n", i, g.Length())
		}
	}

	return ret
}
