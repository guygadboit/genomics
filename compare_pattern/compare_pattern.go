/*
Just pick a pattern out of a big alignment and count how many genomes have what
*/
package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"os"
	"log"
)

func countPattern(g *genomes.Genomes, start, end int) map[string][]string {
	ret := make(map[string][]string)
	var env genomes.Environment

	err := env.Init(g, start, end-start, 0)
	if err != nil {
		log.Fatal(err)
	}

	for i := 0; i < g.NumGenomes(); i++ {
		pat := g.Nts[i][start:end]

		silent, differences := env.Replace(pat)
		if !silent {
			continue
		}

		spat := fmt.Sprintf("%s (%d)", string(pat), differences)
		record, _ := ret[spat]
		ret[spat] = append(record, g.Names[i])
	}
	return ret
}

func main() {
	var start, end int
	var orfs string

	// 22922:22927 is the range you are interested in
	flag.IntVar(&start, "s", 22922, "Start (1-based, inclusive)")
	flag.IntVar(&end, "e", 22927, "End (1-based, inclusive)")
	flag.StringVar(&orfs, "orfs", "", "ORFS file")
	flag.Parse()

	argi := len(os.Args) - flag.NArg()
	g := genomes.LoadGenomes(os.Args[argi], orfs, false)
	g.RemoveGaps()

	start -= 1	// because 1-based
	results := countPattern(g, start, end)

	for k, v := range results {
		fmt.Printf("%s is in the following:\n", k)
		for _, name := range v {
			fmt.Println(name)
		}
	}
}
