package main

import (
	"genomics/genomes"
	"fmt"
)

func main() {
	g := genomes.LoadGenomes("../fasta/more_relatives.fasta",
		"../fasta/WH1.orfs", false)
	g.RemoveGaps()

	/*
	translation := genomes.Translate(g, 0)

	for pos := 0; pos != -1; {
		pos = translation.Find([]byte("DSKVG"), pos)
		if pos == -1 {
			break
		}
		fmt.Printf("%d\n", pos)
		pos++
	}
	*/

	for i := 0; i < g.NumGenomes(); i++ {
		var env genomes.Environment
		err := env.Init(g, 22885, 100, i)
		if err != nil && err.Error() == "Not in ORF" {
			continue
		}
		fmt.Printf(">%s\n", g.Names[i])
		fmt.Printf("%s\n", string(env.ProteinShort()))
	}
}
