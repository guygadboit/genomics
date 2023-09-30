package main

import (
	"fmt"
	"genomics/genomes"
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
	// g := genomes.LoadGenomes("/fs/f/genomes/human/human.fasta.gz", "", true)
	fmt.Printf("Loaded\n")
	index.Build(g, "./index", 3, true)
	// index.Build(g, "/fs/f/genomes/human/index", 6, true)
	index.Save()
}

func testUseIndex() {
	// pattern := []byte("TTTTTTTC")
	pattern := []byte("ATCT")
	// pattern := []byte("CTCCTCGGCGG")

	var s genomes.IndexSearch
	for s.Init("./index", pattern); !s.End(); s.Next() {
	// for s.Init("/fs/f/genomes/human/index", pattern); !s.End(); s.Next() {
		pos, _ := s.Get()
		fmt.Printf("%d\n", pos)
	}

	// Check you can run it twice
	/*
	for s.Start(); !s.End(); s.Next() {
		pos, _ := s.Get()
		fmt.Printf("restarted: %d\n", pos)
	}
	*/
}

func main() {
	// testBuildIndex()
	// testProfile()
	testUseIndex()
}
