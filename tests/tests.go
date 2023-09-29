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
	g := genomes.LoadGenomes("../fasta/WH1.fasta", "", false)
	// g := genomes.LoadGenomes("/fs/f/genomes/human/human.fasta.gz", "", false)
	fmt.Printf("Loaded\n")
	index.Build(g, "./index", 6, true)
	index.Save()
}

func testUseIndex() {
	g := genomes.LoadGenomes("../fasta/WH1.fasta", "", false)
	var s genomes.IndexSearch
	for s.Init("./index", 6, g, 0, []byte("CTCCTCGGCGGG")); !s.End(); s.Next() {
		fmt.Println(s.Get())
	}
}

func main() {
	// testBuildIndex()
	// testProfile()
	testUseIndex()
}
