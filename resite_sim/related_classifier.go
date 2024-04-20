package main

import (
	"fmt"
	"genomics/genomes"
	"reflect"
)

type Classifier struct {
	relatives *genomes.Genomes
}

/*
Load some bunch of relatives to compare to. This works pretty well if you just
use the BANALs
*/
func (c *Classifier) Init() {
	c.relatives = genomes.LoadGenomes("../fasta/CloseRelatives.fasta",
		"../fasta/WH1.orfs", false)
	/*
	c.relatives = genomes.LoadGenomes("../fasta/ACCRealigned.fasta",
		"../fasta/WH1.orfs", false)
	*/
	c.relatives.RemoveGaps()
}

/*
g should contain one genome that shares the alignments of the Classifier's
relatives. Return the number of RE sites that match at least one relative.
exclude is the relative to exclude (because the mutant itself is a mutation of
that genome).
*/
func (c *Classifier) CountShared(g *genomes.Genomes, exclude int) int {
	_, _, _, _, sites := FindRestrictionMap(g)

	matched := 0
	for _, site := range sites {
		for i := 0; i < c.relatives.NumGenomes(); i++ {
			if i == exclude {
				continue
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
Create a mutant based on reference and see if the classifier can tell that it's
a mutant. Do this numTrials times and return the number of successes
*/
func (c *Classifier) OneTrial(reference int, numTrials int) int {
	success := 0

	/*
	First see how many sites are shared between the one we're going to "tamper"
	with and the relatives. If fewer than this number are shared afterwards
	then we know tampering has happened.
	*/
	mutant := c.relatives.Filter(reference)
	benchmark := c.CountShared(mutant, reference)

	for i := 0; i < numTrials; i++ {
		mutant := c.relatives.Filter(reference)
		mutant.DeepCopy(0)
		Tamper(mutant, RE_SITES, 1, 1)

		count := c.CountShared(mutant, reference)
		if count < benchmark {
			success++
			fmt.Printf("Success %d/%d\n", success, i+1)
		}
	}
	fmt.Printf("Success rate %d/%d\n", success, numTrials)
	return success
}

func (c *Classifier) Trial(numTrials int) int {
	success := 0
	n := c.relatives.NumGenomes()
	for i := 0; i < n; i++ {
		success += c.OneTrial(i, numTrials)
	}
	fmt.Printf("Overall success rate %d/%d (%.2f%%)\n", success, numTrials*n,
		float64(success*100)/float64(numTrials*n))
	return success
}
