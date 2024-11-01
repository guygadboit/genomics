package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
)

func main() {
	flag.Parse()

	for _, arg := range flag.Args() {
		g := genomes.LoadGenomes(arg, "", true)
		fmt.Println(g.Names[0])
		var alphabet string
		if g.IsProtein() {
			alphabet = mutations.AA_ALPHABET
		} else {
			alphabet = mutations.NT_ALPHABET
		}
		nd := mutations.NewNucDistro(mutations.NewGenomeIterator(g), alphabet)
		nd.Show()
	}
}
