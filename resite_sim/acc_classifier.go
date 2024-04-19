package main

import (
	"fmt"
	"genomics/genomes"
	"reflect"
)

type Classifier struct {
	relatives *genomes.Genomes
}

func (c *Classifier) Init() {
	c.relatives = genomes.LoadGenomes("../fasta/ACC_Alignment.fasta",
		"../fasta/WH1.orfs", false)
	c.relatives.RemoveGaps()
}

func (c *Classifier) IsTampered(g *genomes.Genomes) bool {
	_, _, _, _, sites := FindRestrictionMap(g)

	// Look for whether we can find the sites in any of the relatives. The
	// problem here is that g and relatives aren't aligned. You might 
	matched := 0
	for _, site := range sites {
		for i := 0; i < c.relatives.NumGenomes(); i++ {
			if reflect.DeepEqual(c.relatives.Nts[i][site:site+6],
				g.Nts[0][site:site+6]) {
				matched++
				break
			}
		}
	}
	return matched > 5
}

func (c *Classifier) Trial(numTrials int) int {
	success := 0
	for i := 0; i < numTrials; i++ {
		mutant := c.relatives.Filter(0)
		mutant.DeepCopy(0)
		Tamper(mutant, RE_SITES, 3, 3)
		looksTampered := c.IsTampered(mutant)
		if looksTampered {
			fmt.Println("Success")
			success++
		} else {
			fmt.Println("Fail")
		}
	}
	fmt.Printf("Success rate %d/%d\n", success, numTrials)
	return success
}

var classifier Classifier

func GetClassifier() *Classifier {
	return &classifier
}

func init() {
	classifier.Init()
}
