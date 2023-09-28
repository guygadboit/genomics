package main

import (
	"genomics/genomes"
)

func testProfile() {
	s := "TCATGACGTTCGTGTTGTTTTAATCTAAACG"
	prof := genomes.CalcProfile([]byte(s))
	prof.Show()
}

func testIndex() {
	g := LoadGenomes("../fasta/WH1.fasta", "", false)
}

func main() {
	testProfile()
}
