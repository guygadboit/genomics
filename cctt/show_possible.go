package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
)

func ShowPossible() {
	g := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	muts := mutations.PossibleSilentMuts(g, 0)
	highlights := make([]genomes.Highlight, 0)
	for _, mut := range muts {
		var char byte
		if mut.From == 'T' && mut.To == 'C' {
			char = 'V'
		} else {
			continue
			char = 'v'
		}
		highlights = append(highlights,
			genomes.Highlight{mut.Pos, mut.Pos + 1, char})
	}
	g.SaveClu("possible-silent.clu", highlights, 0)
	fmt.Printf("%d locations marked. Wrote possible-silent.clu\n",
		len(highlights))
}
