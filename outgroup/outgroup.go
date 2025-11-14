package outgroup

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
)

// Return the number of differences. Assume a and b are the same length
func CompareRegion(a []byte, b []byte, start, end int) int {
	aRegion := a[start:end]
	bRegion := b[start:end]

	ret := 0
	for i, val := range aRegion {
		if bRegion[i] != val {
			ret++
		}
	}

	return ret
}

// Describes proximity to a relative
type Proximity struct {
	Which       int // Which relative
	Differences int // How many differences
	Window      int // At what window size
}

func ShowProximities(g *genomes.Genomes, which int,
	sitePos int, proximities []Proximity) {
	for _, p := range proximities {
		fmt.Println(g.Names[p.Which], p.Window, p.Differences)
	}
}

type Set map[int]bool

// Return the relatives that are closest window either side of the site
func FindClosest(g *genomes.Genomes, which int,
	sitePos int, siteSize int, window int, exclude Set) []Proximity {
	ret := make([]Proximity, 0)

	for i := 0; i < g.NumGenomes(); i++ {
		if i == which {
			continue
		}
		if exclude[i] {
			continue
		}
		end := sitePos
		start := max(0, end-window)
		differences := CompareRegion(g.Nts[which], g.Nts[i], start, end)

		start = min(g.Length(), sitePos+siteSize)
		end = min(g.Length(), start+window)

		differences += CompareRegion(g.Nts[which], g.Nts[i], start, end)
		ret = append(ret, Proximity{i, differences, window})
	}
	utils.SortByKey(ret, true, func(p Proximity) int {
		return p.Differences
	})
	return ret
}

type Algo int

// KEEP_BEST is better. ORIGINAL is what Hensel used in his paper.
const (
	ORIGINAL = iota
	KEEP_BEST
)

// Find the num unambiguous closest either side of the site
func FindNumClosest(g *genomes.Genomes, which int,
	sitePos int, siteSize int, window int, num int, algo Algo) []Proximity {
	num = min(num, g.NumGenomes()-1)
	ret := make([]Proximity, 0, num)
	got := make(Set)

	for {
		need := num - len(ret) // how many more relatives do we need?

		// Find the closest ones, excluding those we already have
		prox := FindClosest(g, which, sitePos, siteSize, window, got)

		// Take any clear leaders
		for _, p := range prox {
			if p.Differences < prox[need].Differences {
				ret = append(ret, p)
				got[p.Which] = true
			}
		}

		if len(ret) == num {
			return ret
		}

		if algo == ORIGINAL {
			/*
				Under the original algorithm, if we didn't find the number we
				wanted (which was three), we threw them all away and then
				widened the window. This meant you threw away some quite good
				matches and ended up with worse ones over a longer window.
			*/
			ret = ret[:0]
			got = make(Set)
		}
		window++
	}
}
