package main

import (
	"fmt"
	"genomics/genomes"
)

/*
	I wonder if any of the inserts appear in the genomes inside larger repeats?
	Let's find out. For this we use an index and the actual genome (so we can
	backtrack either side)
*/
func findSatellites(genome *genomes.Genomes, index string,
	ins *Insertion, name string, verbose bool) int {
	var search genomes.BidiIndexSearch
	nts := genome.Nts[0]
	n := len(ins.nts)
	var ret int

	clamp := func(x int) int {
		if x < 0 {
			return 0
		}
		if x > len(nts) {
			return len(nts)
		}
		return x
	}

outer:
	for search.Init(index, ins.nts); !search.End(); search.Next() {
		pos, _ := search.Get()

		for start := clamp(pos - 3); start < pos; start++ {
			for end := pos + n + 1; end <= clamp(pos+n+3); end++ {
				var search2 genomes.IndexSearch
				for search2.Init(index,
					nts[start:end]); !search2.End(); search2.Next() {
					pos2, _ := search2.Get()
					if pos2 == start {
						continue
					}
					if verbose {
						fmt.Printf("Found %d %s at %d (-%d +%d) "+
							"inside repeat in %s\n",
							ins.id, string(ins.nts), start,
							pos-start, end-(pos+n), name)
					}
					ret += 1
					continue outer
				}
			}
		}
	}
	return ret
}
