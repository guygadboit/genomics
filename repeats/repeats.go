/*
	The idea was to find all the repeats in a given genome. But it's going to
	take ages.
*/
package main

import (
	"fmt"
	"genomics/genomes"
)

type Repeat struct {
	count  int
	length int
}

/*
	Call cb for all the subsequences of a given length. If it returns false
	stop the iteration.
*/
func findSubsequences(genome *genomes.Genomes,
	length int, cb func(int, []byte) bool) {
	nts := genome.Nts[0]
	for i := 0; i < genome.Length() - length; i++ {
		if !cb(i, nts[i:i+length]) {
			break
		}
	}
}

func FindRepeats(genome *genomes.Genomes,
	index string, minLength int, maxLength int) []Repeat {

	ret := make([]Repeat, 0)
	// Tracks which positions are in repeats
	inRepeat := make([]bool, genome.Length())

	markInRepeat := func(start, length int) {
		for i := start; i < start + length; i++ {
			inRepeat[i] = true
		}
	}

	for l := maxLength; l >= minLength; l-- {
		fmt.Printf("Searching at length %d\n", l)
		rep := Repeat{length: l}
		findSubsequences(genome, l, func(pos int, nts []byte) bool {
			var search genomes.IndexSearch
			for search.Init(index, nts); !search.End(); search.Next() {
				repeatPos, _ := search.Get()
				if repeatPos == pos || inRepeat[repeatPos] {
					continue
				}
				fmt.Printf("%s\n", string(nts[pos:pos+l]))
				rep.count++
				markInRepeat(repeatPos, l)
			}
			if rep.count > 0 {
				markInRepeat(pos, l)
			}
			return true
		});
		if rep.count > 0 {
			ret = append(ret, rep)
		}
	}
	return ret
}

func main() {
	g := genomes.LoadGenomes("/fs/f/genomes/bacteria/"+
		"Legionella/Legionella.fasta", "", true)
	index := "/fs/f/genomes/bacteria/Legionella/index"
	repeats := FindRepeats(g, index, 12, 200)
	for i := 0; i < len(repeats); i++ {
		fmt.Printf("%d, %d\n", repeats[i].count, repeats[i].length)
	}
}
