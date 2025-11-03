package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"genomics/utils"
	"math/rand"
)

func CalcSimilarities(g *genomes.Genomes, which int) []float64 {
	ret := make([]float64, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret[i] = g.SequenceSimilarity(which, i, false)
	}
	return ret
}

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

type Proximity struct {
	which       int
	differences int
	window      int
}

// We assume Proximity is sorted by fewest differences first
func HaveUniqueThreeBest(p []Proximity) bool {
	// We just need fresh air between the third and the fourth to satisfy the
	// criterion.
	return p[3].differences > p[2].differences
}

func FindClosest(g *genomes.Genomes, which int,
	sitePos int, siteSize int, window int) []Proximity {
	ret := make([]Proximity, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		if i == which {
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
		return p.differences
	})
	if !HaveUniqueThreeBest(ret) {
		return FindClosest(g, which, sitePos, siteSize, window+1)
	}
	/*
		for _, p := range ret {
			fmt.Println(sitePos, window, g.Names[p.which], p)
		}
	*/
	return ret
}

func CompareRelatives(g *genomes.Genomes, which int, ss []float64) {
	for _, site := range hotspots.RE_SITES {
		for search := genomes.NewLinearSearch(g,
			0, site, 0.0); !search.End(); search.Next() {
			pos, _ := search.Get()
			closest := FindClosest(g, which, pos, len(site), 50)
			// Now the question is in how many of the three closest does the
			// actual site alo match?
			count := 0
			for _, p := range closest[:3] {
				differences := CompareRegion(g.Nts[which], g.Nts[p.which],
					pos, pos+len(site))
				if differences == 0 {
					count++
				}
				fmt.Printf("%s at %d in %s (closest over %dnts either "+
					"side) has %d differences in the site\n",
					string(site), pos, g.Names[p.which], p.window, differences)
			}
			fmt.Printf("%d/3 are completely the same\n", count)
		}
	}
}

func AddSite(g *genomes.Genomes, pattern []byte, which int) {
	for {
		pos := rand.Intn(g.Length() - len(pattern))
		for i := pos; i < g.Length()-len(pattern); i++ {
			alreadyThere := true
			for j := 0; j < len(pattern); j++ {
				if g.Nts[which][i+j] != pattern[j] {
					alreadyThere = false
				}
			}
			if alreadyThere {
				continue
			}
			silent, numMuts, _ := genomes.IsSilentWithReplacement(g, i,
				0, 0, pattern)
			if silent && numMuts == 1 {
				fmt.Printf("New %s site at %d\n", string(pattern), i)
				copy(g.Nts[which][i:i+len(pattern)], pattern)
				// g.Save("Modified", "modified.fasta", 0)
				return
			}
		}
	}
}

func main() {
	var tamper bool

	flag.BoolVar(&tamper, "tamper", false, "Whether to adjust sites")
	flag.Parse()

	g := genomes.LoadGenomes("../fasta/Hassanin.fasta",
		"../fasta/WH1.orfs", false)

	if tamper {
		for _, site := range hotspots.RE_SITES {
			AddSite(g, site, 0)
		}
	}

	ss := CalcSimilarities(g, 0)
	CompareRelatives(g, 0, ss)
}
