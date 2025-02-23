package main

import (
	"fmt"
	"flag"
	"genomics/genomes"
	"genomics/utils"
)

func main() {
	var outname string

	flag.StringVar(&outname, "o", "dealigned.fasta", "output name")
	flag.Parse()

	fname := flag.Args()[0]
	g := genomes.LoadGenomes(fname, "", false)
	dealigned := g.Dealign()

	fd, fp := utils.WriteFile(outname)
	defer fd.Close()

	for _, g := range dealigned {
		fmt.Fprintf(fp, ">%s\n", g.Names[0])
		utils.Wrap(fp, g.Nts[0])
	}
	fp.Flush()

	fmt.Printf("Wrote %s\n", outname)
}
