package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"math/rand"
	"os"
	"path/filepath"
	"strings"
)

// The alleles in a given location
type Alleles struct {
	Pos int
	Nts map[string][]int // For the nts in a given codon, which genomes have it?
	Aas map[byte][]int   // For a given AA, which genomes have it?
}

func NewAlleles(pos int) *Alleles {
	var ret Alleles
	ret.Pos = pos
	ret.Nts = make(map[string][]int)
	ret.Aas = make(map[byte][]int)
	return &ret
}

func (a *Alleles) Add(codon *genomes.Codon, genomeIndex int) {
	_, there := a.Nts[codon.Nts]
	if !there {
		a.Nts[codon.Nts] = make([]int, 0)
	}
	a.Nts[codon.Nts] = append(a.Nts[codon.Nts], genomeIndex)

	_, there = a.Aas[codon.Nts]
	if !there {
		a.Aas[codon.Aa] = make([]int, 0)
	}
	a.Aas[codon.Aa] = append(a.Aas[codon.Aa], genomeIndex)
}

func Translate(g *genomes.Genomes) []genomes.Translation {
	ret := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret[i] = genomes.Translate(g, i)
	}
	return
}

func GetAlleles(translations []genomes.Translation, codonPos int) *Alleles {
	ret := NewAlleles(codonPos)
	for i, _ := range translations {
		codon := translations[i][codonPos]
		ret.Add(&codon)
	}
	return ret
}

type Result struct {
	UniqueAas                  int
	SoleOutlierAas             int
	UniqueSilentCodons         int
	UniqueNonSilentCodons      int
	SoleOutlierSilentCodons    int
	SoleOutlierNonSilentCodons int
}

// One fo each genome
type Results []Result

// If there is a unique Aa in there, return the genome that has it. Otherwise
// -1
func (a *Alleles) FindUniqueAa() int {
	for _, v := range a.Aas {
		if len(v) == 1 {
			return v[0]
		}
	}
	return -1
}

func (a *Alleles) FindSoleOutlierAa() int {
	others := make(map[byte]bool)
	ret := -1

	for k, v := range a.Aas {
		if len(v) == 1 {
			ret = v[0]
		} else {
			others[k] = true
			if len(others) > 1 {
				return -1
			}
		}
	}

	return ret
}

func (a *Alleles) FindUniqueNts() int {
	for _, v := range a.Nts {
		if len(v) == 1 {
			return v[0]
		}
	}
	return -1
}

func (a *Alleles) FindSoleOutlierNts() int {
	others := make(map[string]bool)
	ret := -1

	for k, v := range a.Nts {
		if len(v) == 1 {
			ret = v[0]
		} else {
			others[k] = true
			if len(others) > 1 {
				return -1
			}
		}
	}

	return ret
}

type FindFunc func(a *Alleles) int

// Update results with the counts based on the alleles at a particular location
func (a *Alleles) Count(g *genomes.Genomes, results Results) {

	index := a.FindUniqueAa()
	if index != -1 {
		results[index].UniqueAas++
	}

	index = a.FindSoleOutlierAa()
	if index != -1 {
		results[index].SoleOutlierAas++
	}

	index = a.FindUniqueNts()
	// Decide if silent or not... YOU ARE HERE




}
