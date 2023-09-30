package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
)

func BuildIndex(fname string, dir string, length int) {
	fmt.Printf("Loading %s...\n", fname)
	g := genomes.LoadGenomes(fname, "", true)
	fmt.Printf("Loaded\n")

	var index genomes.Index
	index.Build(g, dir, length, true)
	index.Save()
}

func main() {
	var fname, dir string
	var length int

	flag.StringVar(&fname, "file", "", "Fasta file name")
	flag.StringVar(&dir, "dir", "", "Where to put the index")
	flag.IntVar(&length, "length", 6, "Index word length")

	flag.Parse()

	BuildIndex(fname, dir, length)
}
