package main

import (
	"genomics/genomes"
	"sort"
)

/*
	How many nts are the same? Return the highest number that are and what
	value that corresponds to. But ignore anything with a negative value
	because that doesn't count.
*/
func countSame(nts []int) (int, int) {
	// Make a copy of nts with anything < 0 filtered out and then sort it.
	sorted := make([]int, 0, len(nts))
	for i := 0; i < len(nts); i++ {
		if nts[i] >= 0 {
			sorted = append(sorted, nts[i])
		}
	}
	sort.Ints(sorted)

	var count, best, bestVal int

	/*
	 Looping one past the end so that if we have a series of matches right
	 at the end of the list we can handle what happens when that series
	 ends in the same way (the else block in the loop). In other words, we
	 go into the else when we find something different, or we've gone past
	 the end.
	*/
	for i := 1; i < len(sorted)+1; i++ {
		if i < len(sorted) && sorted[i] == sorted[i-1] {
			count++
		} else {
			if count > best {
				best = count
				bestVal = sorted[i-1]
			}
			count = 0
		}
	}

	// +1 because we're counting how many subsequent nts match the first one.
	// So the total number of matches is that +1 (to include that first one).
	return best + 1, bestVal
}

/*
	Do nts appear more than once in genome, and if they do, how long is the
	total match? Return the length of the repeating section and how many times
	it repeats. The length is the longest we find and the count how many
	repeats of any length greater than the pattern length (don't set too much
	store by the number of repeats, the length is more interesting).
*/
func findSatellites(genome *genomes.Genomes, index string,
	pattern []byte, name string) (int, int) {
	nts := genome.Nts[0]
	var search genomes.IndexSearch
	search.Init(index, pattern)
	positions := genomes.SearchAll(&search)

	// Extend the positions as far as we can in the given direction (1 or -1).
	// Return how far and how many.
	extend := func(direction int) (int, int) {
		var distance, count int
		// Using ints because they're much easier to sort in Go
		extensions := make([]int, len(positions))
		for i := 0; i < len(extensions); i++ {
			extensions[i] = -1 // means invalid
		}

		for distance = 1; ; distance++ {
			for i := 0; i < len(extensions); i++ {
				if extensions[i] == -2 {
					continue
				}
				var pos int
				switch direction {
				case 1:
					pos = positions[i] + len(pattern) + distance
				case -1:
					pos = positions[i] - distance
				}
				if pos >= 0 && pos < genome.Length() {
					extensions[i] = int(nts[pos])
				}
			}

			same, val := countSame(extensions)
			//fmt.Println(extensions, same, val)
			if same == 1 {
				break
			}
			if count == 0 {
				count = same
			}

			// OK now any of these positions that don't contain the "consensus"
			// are "bust", which we indicate with -2
			for j := 0; j < len(extensions); j++ {
				if extensions[j] != val {
					extensions[j] = -2
				}
			}
		}
		return distance - 1, count
	}

	prefixLen, prefixCount := extend(-1)
	suffixLen, suffixCount := extend(1)
	return prefixLen + len(pattern) + suffixLen, prefixCount + suffixCount
}
