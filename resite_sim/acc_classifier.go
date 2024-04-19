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
}

func (c *Classifier) IsTampered(g *genomes.Genomes) bool {
	_, _, _, _, sites := FindRestrictionMap(g)

	// Look for whether we can find the sites in any of the relatives. The
	// problem here is that g and relatives aren't aligned. FIXME YOU ARE HERE.
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

	fmt.Println(matched)
	return true
}

var classifier Classifier

func GetClassifier() *Classifier {
	return &classifier
}

func init() {
	classifier.Init()
}
