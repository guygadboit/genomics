package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

func FindTandemRepeats(genome *genomes.Genomes,
	length int, verbose bool) []int {
	ret := make([]int, 0)
	nts := genome.Nts[0]
searching:
	for i := 0; i < genome.Length()-length*2; i++ {
		pattern := nts[i : i+length]
		if utils.IsSilly(pattern) || !utils.IsRegularPattern(pattern) {
			continue
		}
		for j := 0; j < length; j++ {
			if nts[i+j] != nts[i+length+j] {
				continue searching
			}
		}
		ret = append(ret, i)
		if verbose {
			fmt.Printf("%d (%d): %s %s\n", i, length,
				string(nts[i:i+length]), genome.Names[0])
		}
		i += length
	}
	return ret
}

func main() {
	var dealign, verbose bool

	flag.BoolVar(&dealign, "d", false, "Dealign, i.e. assume input is an"+
		"alignment, not something like human which needs merging")
	flag.BoolVar(&verbose, "v", false, "Verbose")
	flag.Parse()

	for _, arg := range flag.Args() {
		var gs []*genomes.Genomes
		input := genomes.LoadGenomes(arg, "", !dealign)
		if dealign {
			gs = input.Dealign()
		} else {
			gs = []*genomes.Genomes{input}
		}

		for _, g := range gs {
			repeats := FindTandemRepeats(g, 3, verbose)
			fmt.Println(g.Names[0], len(repeats))
			/*
				for i := 10; i < 60; i++ {
					FindTandemRepeats(g, i, true)
				}
			*/
		}
	}
}
