package main

import (
	"fmt"
	"genomics/genomes"
	"sort"
)

// How many nts are the same? Return the highest number that are and what value
// that corresponds to. But ignore anything with a negative value because that
// doesn't count.
func countSame(nts []int) (int, int) {
	// Make a copy of nts with anything < 0 filtered out and then sort it.
	sorted := make([]int, 0, len(nts))
	for i := 0; i < len(nts); i++ {
		if nts[i] >= 0 {
			sorted = append(sorted, nts[i])
		}
	}
	sort.Ints(sorted)

	var count, best, bestVal, lastVal int

	// Looping one past the end so that if we have a series of matches right at
	// the end of the list we can handle what happens when that series ends in
	// the same way (the else block in the loop). In other words, we go into
	// the else when we find something different, or we've gone past the end.
	for i := 1; i < len(sorted) + 1; i++ {
		if i < len(sorted) && sorted[i] == sorted[i-1] {
			count += 1
		} else {
			if count > best {
				best = count
				bestVal = sorted[i-1]
			}
			count = 0
		}
	}

	return best, bestVal
}


/*
	Do nts appear more than once in genome, and if they do, how long is the
	total match? Return the length of the repeating section and how many times
	it repeats. The length can be the longest we find and the count how many
	repeats.
*/
func findRepeats(genome *genomes.Genomes, index string,
	pattern []byte, name string, verbose bool) (int, int) {
	nts := genome.Nts[0]
	var search genomes.BidiIndexSearch
	search.Init(index, pattern)
	positions := genomes.SearchAll(&search)

	// Using ints because they're much easier to sort in Go
	prefixes := make([]int, len(positions))
	for i := 0; i < len(prefixes); i++ {
		prefixes[i] = -1	// means invalid
	}

	var distance int
	for distance = 1; ; distance++ {
		for i := 0; i < len(prefixes); i++ {
			pos := positions[i] - distance
			if pos >= 0 {
				if prefixes[i] != 2 {
					prefixes[i] = int(nts[pos])
				}
			}
			same, val := countSame(prefixes)
			if same == 0 {
				break
			}

			// OK now any of these positions that don't contain the "consensus"
			// are "bust", which we indicate with -2
			for j := 0; j < len(prefixes); j++ {
				if prefixes[j] != val {
					prefixes[j] = -2
				}
			}
		}
	}
}
