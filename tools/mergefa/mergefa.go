package main

import (
	"fmt"
	"genomics/genomes"
	"flag"
)

func main() {
	var name, outName string

	flag.StringVar(&name, "name", "MergedFasta", "Output genome name")
	flag.StringVar(&outName, "o", "output.fasta", "Output filename")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", true)
	g.Save(name, outName, 0)
	fmt.Printf("Wrote %s\n", outName)
}
