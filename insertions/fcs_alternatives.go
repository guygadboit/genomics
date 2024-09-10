package main

import (
	"fmt"
	"genomics/genomes"
)

func CheckAlternatives() {
	wh1 := genomes.LoadGenomes("../fasta/WH1.fasta",
		"../fasta/WH1.orfs", false)

	// Now we need to delete the actual FCS because we're going to be putting
	// it back in again with the insertions. If you don't delete it you will
	// get lots of false homology.
	wh1d := wh1.Filter(0)
	nts := make([]byte, wh1.Length()-12)
	copy(nts, wh1.Nts[0][:23600])
	copy(nts[23600:], wh1.Nts[0][23612:])
	wh1d.Nts[0] = nts

	var env genomes.Environment

	// The alternatives have to be found in the real WH1
	env.Init(wh1, 23600, 12, 0)
	alternatives := env.FindAlternatives(12, true)
	fmt.Printf("There are %d alternatives\n", len(alternatives))

	sources := GetSources(BACTERIA)
	insertions := make([]Insertion, 0)

	for i, alt := range alternatives {
		ins := Insertion{i+1, 23601,
			alt.Nts, 2, false, false,
			0, UNKNOWN, false, 50, 50, nil,
		}
		insertions = append(insertions, ins)
	}

	// And let's throw the actual FCS in the list as well
	insertions = appendFCS(insertions)
	data := InsertionData{nil, nil, insertions, nil}

	CountInGenomes(wh1d, &data, sources, nil, 0.0, 1, APPEND)
	OutputMatches(insertions, "alternative-fcs-matches.txt")
}
