package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

func main() {
	var include, exclude, outName string
	var subSeq string

	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.StringVar(&exclude, "e", "", "Genomes to exclude")
	flag.StringVar(&outName, "o", "output.clu", "Output filename")
	flag.StringVar(&subSeq, "r", "", "range")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", false)

	clamp := func(x int) int {
		if x < 0 {
			return 0
		}
		if x > g.Length() {
			return g.Length()
		}
		return x
	}

	if subSeq != "" {
		limits := utils.ParseInts(subSeq, ":")
		start := clamp(limits[0])
		end := clamp(limits[1])
		for i, nts := range g.Nts {
			g.Nts[i] = nts[start:end]
		}
	}

	var which []int
	if include != "" {
		inc := utils.ToSet(utils.ParseInts(include, ","))
		for i := 0; i < g.NumGenomes(); i++ {
			if inc[i] {
				which = append(which, i)
			}
		}
	}

	g.SaveClu(outName, nil, which...)
	fmt.Printf("Wrote %s\n", outName)
}