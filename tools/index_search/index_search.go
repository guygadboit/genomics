package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"log"
	"strings"
)

func Search(root string, needle []byte, g *genomes.Genomes, context int) {
	count := 0
	for search := genomes.NewBidiIndexSearch(root,
		needle); !search.End(); search.Next() {
		pos, err := search.Get()
		if err != nil {
			log.Fatal(err)
		}
		var f string
		if search.IsForwards() {
			f = "+"
		}
		fmt.Printf("%d%s\n", pos+1, f)

		if context != 0 {
			start := max(0, pos-context)
			end := min(g.Length(), pos+len(needle)+context)
			fmt.Printf("%s\n", g.Nts[0][start:end])
		}

		count++
	}
	fmt.Printf("%d matches\n", count)
}

func main() {
	var (
		pattern string
		fasta   string
		context int
	)

	flag.StringVar(&pattern, "p", "", "Pattern")
	flag.IntVar(&context, "c", 0, "Show context around matches")
	flag.StringVar(&fasta, "fasta", "", "Fasta file (needed if context != 0)")
	flag.Parse()

	var g *genomes.Genomes
	if context != 0 {
		g = genomes.LoadGenomes(fasta, "", true)
	}

	pattern = strings.ToUpper(pattern)
	needle := []byte(pattern)
	for _, arg := range flag.Args() {
		fmt.Println(arg)
		Search(arg, needle, g, context)
	}
}
