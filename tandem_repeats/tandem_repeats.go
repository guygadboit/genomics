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

func count(repeats []int) (int, int) {
	roiStart, roiEnd := 3000, 4500
	var roi, others int

	for _, r := range repeats {
		if r >= roiStart && r < roiEnd {
			roi++
		} else {
			others++
		}
	}
	return roi, others
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
			roi, others := count(repeats)
			ratio := float64(roi * (30000/1500))/float64(others)
			fmt.Printf("%s %d %d %.4f\n", g.Names[0], roi, others, ratio)
			/*
				for i := 10; i < 60; i++ {
					FindTandemRepeats(g, i, true)
				}
			*/
		}
	}
}
