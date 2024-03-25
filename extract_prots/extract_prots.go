package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
)

type Range struct {
	start, end int
}

func extractRanges(g *genomes.Genomes, ranges []Range) [][]byte {
	ret := make([][]byte, len(g.Nts))
	for i := 0; i < g.NumGenomes(); i++ {
		for _, r := range ranges {
			ret[i] = append(ret[i],
				genomes.TranslateAlignedShort(g.Nts[i][r.start:r.end])...)
		}
	}
	return ret
}

func main() {
	var fasta string
	var outName string

	flag.StringVar(&fasta, "fasta", "", "fasta file (nucleotides)")
	flag.StringVar(&outName, "o", "output.fasta", "output file")
	flag.Parse()

	ranges := make([]Range, 0)
	for _, arg := range flag.Args() {
		ints := utils.ParseInts(arg, ":")
		if len(ints) != 2 {
			log.Fatalf("Invalid range <%s>", arg)
		}
		// Ranges are input in the usual 1-based convention
		ranges = append(ranges, Range{ints[0] - 1, ints[1]})
	}

	g := genomes.LoadGenomes(fasta, "", false)
	aas := extractRanges(g, ranges)

	var orfs genomes.Orfs
	p := genomes.Genomes{aas, g.Names, orfs}
	p.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
