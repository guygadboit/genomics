package main

import (
	"genomics/genomes"
)

func testProfile() {
	s := "TCATGACGTTCGTGTTGTTTTAATCTAAACG"
	prof := genomes.CalcProfile([]byte(s))
	prof.Show()
}

func main() {
	testProfile()
}
