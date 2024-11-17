package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"flag"
)

func main() {
	flag.Parse()
	fname := flag.Arg(0)
	g := genomes.LoadGenomes(fname, "", false)
	for i := 0; i < g.NumGenomes(); i++ {
		g.Nts[i] = utils.ReverseComplement(g.Nts[i])
	}

	baseName, _ := utils.SplitExt(fname)
	outName := fmt.Sprintf("%s-reversed.fasta", baseName)
	g.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
