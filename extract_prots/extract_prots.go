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

func extractRanges(g *genomes.Genomes, ranges []Range, reverse bool) [][]byte {
	ret := make([][]byte, len(g.Nts))
	for i := 0; i < g.NumGenomes(); i++ {
		for _, r := range ranges {
			seq := g.Nts[i][r.start:r.end]
			if reverse {
				seq = utils.ReverseComplement(seq)
			}
			ret[i] = append(ret[i], genomes.TranslateAlignedShort(seq)...)
		}
	}
	return ret
}

func checkRanges(g *genomes.Genomes, ranges []Range) {
	for _, r := range ranges {
		if r.start < 0 || r.end >= g.Length() {
			log.Fatalf("Invalid range %d:%d (length is %d)",
				r.start+1, r.end, g.Length())
		}
	}
}

func main() {
	var (
		fasta           string
		outName, format string
		removeGaps      bool
		reverse         bool
	)

	flag.StringVar(&fasta, "fasta", "", "fasta file (nucleotides)")
	flag.StringVar(&outName, "o", "output.fasta", "output file")
	flag.StringVar(&format, "format", "fasta", "output format")
	flag.BoolVar(&removeGaps, "g", false, "remove gaps")
	flag.BoolVar(&reverse, "r", false, "reverse")
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
	if removeGaps {
		g.RemoveGaps()
	}

	checkRanges(g, ranges)
	aas := extractRanges(g, ranges, reverse)

	var orfs genomes.Orfs
	p := genomes.Genomes{aas, g.Names, orfs}

	switch format {
	case "fasta":
		p.SaveMulti(outName)
	case "clu":
		p.SaveClu(outName, nil)
	default:
		log.Fatalf("Invalid format: <%s>", format)
	}

	fmt.Printf("Wrote %s\n", outName)
}
