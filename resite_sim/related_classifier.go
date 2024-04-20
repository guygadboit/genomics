package main

import (
	"fmt"
	"genomics/genomes"
	"math"
	"reflect"
	"sort"
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

	c.relatives = genomes.LoadGenomes("../fasta/CloseRelatives.fasta",
		"../fasta/WH1.orfs", false)

	/*
		c.relatives = genomes.LoadGenomes("../fasta/ACCRealigned.fasta",
		"../fasta/WH1.orfs", false)
	*/
	c.relatives.RemoveGaps()

	sites := make(map[int]bool)
	for i := 0; i < c.relatives.NumGenomes(); i++ {
		_, _, _, _, s := FindRestrictionMap(c.relatives.Filter(i))
		for _, site := range s {
			sites[site] = true
		}
	}
	c.sites = make([]int, 0, len(sites))
	for k, _ := range sites {
		c.sites = append(c.sites, k)
	}
	sort.Ints(c.sites)
	fmt.Println(c.sites)
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
What is the minimum number of shared sites between the set of genomes,
excluding reference?
*/
func (c *Classifier) findBenchmark(reference int) int {
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

	if isTampered(original) {
		ret.FalsePositives++
		if reference == 0 {
			fmt.Println("0th genome reported as tampered")
		}
	} else {
		ret.TrueNegatives++
		if reference == 0 {
			fmt.Println("0th genome reported as legit")
		}
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
