package main

import (
	"genomics/genomes"
)

type NamedOrf struct {
	genomes.Orf
	name string
}

type NamedOrfs []NamedOrf

func (n NamedOrfs) GetName(pos int) string {
	for _, no := range n {
		if pos >= no.Start && pos < no.End {
			return no.name
		}
	}
	return ""
}

var OrfNames NamedOrfs

func init() {
	names := []string{
		"ORF1ab",
		"ORF1ab",
		"S",
		"ORF3a",
		"E",
		"M",
		"ORF6",
		"ORF7a",
		"ORF7b",
		"ORF8",
		"N",
		"ORF10",
	}

	orfs := genomes.LoadOrfs("../fasta/WH1.orfs")
	OrfNames = make(NamedOrfs, len(orfs))

	for i, name := range names {
		OrfNames[i] = NamedOrf{orfs[i], name}
	}
}
