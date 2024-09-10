package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
)

func main() {
	var outName string
	flag.StringVar(&outName, "o", "output.fasta", "Output filename")

	flag.Parse()

	g := genomes.LoadClu(flag.Arg(0), "")
	g.SaveMulti(outName)
	fmt.Printf("Wrote %s\n", outName)
}
