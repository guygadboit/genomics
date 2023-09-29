package main

import (
	"genomics/genomes"
	"fmt"
)

func testProfile() {
	s := "TCATGACGTTCGTGTTGTTTTAATCTAAACG"
	prof := genomes.CalcProfile([]byte(s))
	prof.Show()
}

func testBuildIndex() {
	var index genomes.Index
	fmt.Printf("Loading...\n")
	// g := genomes.LoadGenomes("../fasta/WH1.fasta", "", false)
	g := genomes.LoadGenomes("/fs/f/genomes/human/human.fasta.gz", "", true)
	fmt.Printf("Loaded\n")
	index.Build(g, "/fs/f/genomes/human/index", 6, true)
	index.Save()
}

func testUseIndex() {
	// g := genomes.LoadGenomes("../fasta/WH1.fasta", "", false)
	fmt.Printf("Loading...\n")
	g := genomes.LoadGenomes("/fs/f/genomes/human/human.fasta.gz", "", true)
	fmt.Printf("Loaded\n")

	nts := g.Nts[0]
	pattern := []byte("CTCCTCGGCGGG")

	var s genomes.IndexSearch
	for s.Init("/fs/f/genomes/human/index", 6, g, 0, pattern);
		!s.End(); s.Next() {
		pos, _ := s.Get()
		fmt.Printf("%d: %s\n", pos, string(nts[pos:pos+len(pattern)]))
	}

	/*
	var s2 genomes.Search
	for s2.Init(g, 0, pattern, 0.0); !s2.End(); s2.Next() {
		pos, _ := s2.Get()
		fmt.Printf("check %d: %s\n", pos, string(nts[pos:pos+len(pattern)]))
	}
	*/
}

func main() {
	// testBuildIndex()
	// testProfile()
	testUseIndex()
}
