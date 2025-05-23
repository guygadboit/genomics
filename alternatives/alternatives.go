package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

func main() {
	g := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	pos := 10443

	copy(g.Nts[0][pos:pos+6], []byte("GAGACC"))
	fmt.Printf("Original: %s\n", string(g.Nts[0][pos:pos+6]))

	var nt utils.NtIterator
	for nt.Init(6); !nt.End(); nt.Next() {
		replacement := nt.Get()
		silent, numMuts, err := genomes.IsSilentWithReplacement(g, pos,
			0, 0, replacement)

		if err == nil && silent && numMuts > 0 {
			fmt.Println(numMuts, string(replacement))
		}
	}
}
