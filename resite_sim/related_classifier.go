package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/hotspots"
	"math"
	"reflect"
	"sort"
	"strings"
)

type Classifier struct {
	relatives *genomes.Genomes
	sites     []int // Anywhere there is a site in any of the genomes
}

/*
Load some bunch of relatives to compare to. This works pretty well if you just
use the BANALs
*/
func (c *Classifier) Init() {
	/*
		c.relatives = genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
			"../fasta/WH1.orfs", false)
	*/

	/*
		c.relatives = genomes.LoadGenomes("../fasta/CloseRelatives.fasta",
			"../fasta/WH1.orfs", false)
	*/

	/*
		c.relatives = genomes.LoadGenomes("../fasta/ACCRealigned.fasta",
			"../fasta/WH1.orfs", false)
	*/

	c.relatives = genomes.LoadGenomes("../fasta/Temmam.fasta",
		"../fasta/WH1.orfs", false)

	/*
		c.relatives = genomes.LoadGenomes("../fasta/HassaninPlus.fasta",
			"../fasta/WH1.orfs", false)
	*/

	c.relatives.RemoveGaps()

	// Find the places where any of them have a site
	sites := make(map[int]bool)
	for i := 0; i < c.relatives.NumGenomes(); i++ {
		_, _, _, _, s, _ := FindRestrictionMap(c.relatives.Filter(i))
		for _, site := range s {
			sites[site] = true
		}
	}
	c.sites = make([]int, 0, len(sites))
	for k, _ := range sites {
		c.sites = append(c.sites, k)
	}
	sort.Ints(c.sites)

	/*
		fmt.Println("Places where sites appear in any of the set:")
		fmt.Println(c.sites)
	*/
}

/*
g should contain one genome that shares the alignments of the Classifier's
relatives. Return the number of locations where sites of interest match at
least one relative exclude is the relative to exclude (because the mutant
itself is a mutation of that genome).
*/
func (c *Classifier) CountShared(g *genomes.Genomes, exclude ...int) int {
	matched := 0
	for _, site := range c.sites {
	genomes:
		for i := 0; i < c.relatives.NumGenomes(); i++ {
			for _, x := range exclude {
				if i == x {
					continue genomes
				}
			}
			if reflect.DeepEqual(c.relatives.Nts[i][site:site+6],
				g.Nts[0][site:site+6]) {
				matched++
				break
			}
		}
	}
	return matched
}

/*
What is the actual number of shared sites between reference and the rest?
*/
func (c *Classifier) findBenchmark(reference int) int {
	mutant := c.relatives.Filter(reference)
	return c.CountShared(mutant, reference)
}

/*
What is the minimum number of shared sites between the set of genomes,
excluding reference?
*/
func (c *Classifier) findBenchmark2(reference int) int {
	ret := math.MaxInt
	for i := 0; i < c.relatives.NumGenomes(); i++ {
		if i == reference {
			continue
		}

		mutant := c.relatives.Filter(i)
		// Count how many the others share with i. This means excluding both
		// reference and i
		numShared := c.CountShared(mutant, reference, i)

		if numShared < ret {
			ret = numShared
		}
	}
	return ret
}

type TestResult struct {
	TruePositives  int
	FalsePositives int
	TrueNegatives  int
	FalseNegatives int
}

func (t *TestResult) Combine(other *TestResult) {
	t.TruePositives += other.TruePositives
	t.FalsePositives += other.FalsePositives
	t.TrueNegatives += other.TrueNegatives
	t.FalseNegatives += other.FalseNegatives
}

func (t *TestResult) Show() {
	totalPositives := t.TruePositives + t.FalseNegatives
	totalNegatives := t.TrueNegatives + t.FalsePositives
	sensitivity := float64(t.TruePositives) / float64(totalPositives)
	specificity := float64(t.TrueNegatives) / float64(totalNegatives)

	accuracy := (sensitivity + specificity) / 2

	fmt.Printf("Accuracy %.2f%% Sensitivity: %.2f%% Specificity: %.2f%%\n",
		accuracy*100,
		sensitivity*100,
		specificity*100)
}

/*
Create a mutant based on reference and see if the classifier can tell that it's
a mutant. Do this numTrials times and return the number of successes
*/
func (c *Classifier) OneTrial(reference int, numTrials int) TestResult {
	var ret TestResult
	benchmark := c.findBenchmark(reference)

	isTampered := func(g *genomes.Genomes) bool {
		count := c.CountShared(g, reference)
		return count < benchmark
	}

	original := c.relatives.Filter(reference)

	count := c.CountShared(original, reference)
	fmt.Printf("%s %d/%d shared. Benchmark %d. Tampered: %t\n",
		c.relatives.Names[reference],
		count, len(c.sites), benchmark, count < benchmark)

	if isTampered(original) {
		ret.FalsePositives++
	} else {
		ret.TrueNegatives++
	}

	for i := 0; i < numTrials; i++ {
		mutant := c.relatives.Filter(reference)
		mutant.DeepCopy(0)
		Tamper(mutant, RE_SITES, 3, 3)
		if isTampered(mutant) {
			ret.TruePositives++
		} else {
			ret.FalseNegatives++
		}
	}
	return ret
}

func (c *Classifier) Trial(numTrials int) TestResult {
	var ret TestResult
	n := c.relatives.NumGenomes()
	for i := 0; i < n; i++ {
		r := c.OneTrial(i, numTrials)
		r.Show()
		ret.Combine(&r)
	}
	fmt.Println("Overall result:")
	ret.Show()
	return ret
}

type Matches map[int][]int

func (c *Classifier) ShowMatches(locations []int) Matches {
	ret := make(Matches)
	g := c.relatives
	for _, loc := range locations {
		matched := make([]string, 0)
		ret[loc] = make([]int, 0)
		for i := 1; i < g.NumGenomes(); i++ {
			if reflect.DeepEqual(g.Nts[0][loc:loc+6],
				g.Nts[i][loc:loc+6]) {
				matched = append(matched, fmt.Sprintf("%d", i))
				ret[loc] = append(ret[loc], i)
			}
		}

		fmt.Printf("%5d: WH1 matches %s\n", loc, strings.Join(matched, ","))
	}
	return ret
}

func (c *Classifier) ExploreNeighbours() {
	// These are the 11 places where either SC2 or a close neighbour has a site
	special := []int{2192, 9750, 10443, 11647, 17328,
		17971, 22921, 22922, 23291, 24101, 24508}

	fmt.Println("The 11 special locations:")
	c.ShowMatches(special)

	before := make([]int, len(special))
	for i, v := range special {
		before[i] = v - 6
	}
	fmt.Println("The locations just before them:")
	c.ShowMatches(before)

	after := make([]int, len(special))
	for i, v := range special {
		after[i] = v + 6
	}
	fmt.Println("The locations just after them:")
	c.ShowMatches(after)

	highlights := make([]genomes.Highlight, len(special))
	for i, s := range special {
		highlights[i] = genomes.Highlight{s, s + 6, 'v'}
	}

	c.relatives.SaveClu("alignment.clu", highlights)
	fmt.Println("Wrote alignment.clu")
}

type Site struct {
	Pos   int     // where it is
	Count int     // how many relatives it's in
}

func CountCombinations(sites []Site) {
	n := len(sites)
	numCombinations := 1 << n
	if n > 64 {
		fmt.Println("Not enough bits")
		return
	}

	var (
		total, totalGood int
		totalSites       int
	)

	fmt.Printf("%d loci\n", n)
	for i := 0; i < 1<<n; i++ {
		var numSites, maxGap, weight int

		for j := 0; j < n; j++ {
			if i&(1<<j) != 0 {
				weight += sites[j].Count
				numSites++
				if j > 0 {
					gap := sites[j].Pos - sites[j-1].Pos
					maxGap = max(maxGap, gap)
				}
			}
		}
		total += weight
		if numSites == 6 && maxGap <= 8000 {
			totalGood += weight
		}
		totalSites += numSites
	}
	fmt.Printf("%d/%d %.4f\n", totalGood,
		total, float64(totalGood)/float64(total))

	meanSites := float64(totalSites) / float64(numCombinations)
	fmt.Printf("Mean number of sites: %.02f\n", meanSites)
}

func (c *Classifier) ExploreMissingSites() []Site {
	ret := make([]Site, 0)
	g := c.relatives
	// g = g.FilterOut(0)

	countWithPat := func(location int, pat []byte) int {
		count := 0
		for i := 1; i < g.NumGenomes(); i++ {
			if reflect.DeepEqual(g.Nts[i][location:location+6], pat) {
				count++
			}
		}
		return count
	}

	for _, location := range c.sites {
		have := g.Nts[0][location : location+6]

		// Show places where either SC2 has a site...
		good := false
		var haveSite string
		if hotspots.IsReSite(have) {
			haveSite = "*"
			good = true
		}

		var env genomes.Environment
		env.Init(g, location, 6, 0)
		alternatives := env.FindAlternatives(1, true)

		if !good {
			// ... or where it differs silently from some alternative that has
			// one.
			for _, alt := range alternatives {
				if hotspots.IsReSite(alt.Nts) {
					good = true
					break
				}
			}
			if !good {
				continue
			}
		}

		s := countWithPat(location, have)

		fmt.Printf("Location %d have %s%s shared by %d (%s)\n",
			location+1, string(have), haveSite, s, string(env.Protein()))
		ret = append(ret, Site{location, s})

		for _, alt := range alternatives {
			var isSite string
			if hotspots.IsReSite(alt.Nts) {
				isSite = "*"
				good = true
			}
			count := countWithPat(location, alt.Nts)
			fmt.Printf("%s %d [%d] %s\n",
				string(alt.Nts), count, alt.CodonPos, isSite)
		}
	}
	return ret
}
