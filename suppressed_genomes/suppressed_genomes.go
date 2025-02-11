package main

import (
	"fmt"
	"genomics/genomes"
	"path/filepath"
	"regexp"
)

var ROOT string = "/fs/f/genomes/viruses/suppressed_genomes/"
var DIRS []string = []string{"RdRp", "Spikes", "ORF8"}

type Gene int

const (
	RdRp Gene = iota
	S
	ORF8
)

func GeneString(g Gene) string {
	genes := []string{"RdRP", "S", "ORF8"}
	return genes[g]
}

type Accession struct {
	fname  string
	genome *genomes.Genomes
	gene   Gene

	// The genome name minus the bit that tells you it's ORF8 or whatever. So
	// we can match them up.
	species string
	location string
}

func LoadAll() []Accession {
	pat := regexp.MustCompile(`.*coronavirus strain ([^\s]+).*$`)
	locPat := regexp.MustCompile(`^.*_(.*)`)

	ret := make([]Accession, 0)
	for i, dir := range DIRS {
		fullPath := filepath.Join(ROOT, dir, "MH*.fasta")
		matches, _ := filepath.Glob(fullPath)
		for _, m := range matches {
			g := genomes.LoadGenomes(m, "", false)

			name := g.Names[0]
			matches := pat.FindSubmatch([]byte(name))
			species := string(matches[1])

			matches = locPat.FindSubmatch(matches[1])
			location := string(matches[1])
			fmt.Println(location)

			acc := Accession{m, g, Gene(i), species, location}
			ret = append(ret, acc)
		}
	}
	return ret
}

func CheckLengths(accessions []Accession) {
	for _, acc := range accessions {
		l := acc.genome.Length()
		fmt.Println(acc.gene, acc.genome.Length())
		if l%3 != 0 {
			fmt.Println("NOT Mul3")
		}
	}
	// OK they're all MUL3, So don't need an orfs file. Can just merge them all
	// together.
}

func BySpecies(accessions []Accession) map[string][]Accession {
	ret := make(map[string][]Accession)
	for _, acc := range accessions {
		_, there := ret[acc.species]
		if !there {
			ret[acc.species] = make([]Accession, 0)
		}
		ret[acc.species] = append(ret[acc.species], acc)
	}
	return ret
}

func WhatsMissing(bySpecies map[string][]Accession) {
	for k, v := range bySpecies {
		var gotS, gotO8, gotR bool
		for _, acc := range v {
			switch acc.gene {
			case S:
				gotS = true
			case ORF8:
				gotO8 = true
			case RdRp:
				gotR = true
			}
		}
		if !gotS {
			fmt.Println(k, "missing S")
		}
		if !gotO8 {
			fmt.Println(k, "missing ORF8")
		}
		if !gotR {
			fmt.Println(k, "missing RdRp")
		}
	}
}

func Assemble(bySpecies map[string][]Accession) {
	for k, v := range bySpecies {
		nts := make([]byte, 0)
		for _, acc := range v {
			nts = append(nts, acc.genome.Nts[0]...)
		}

		g := genomes.NewGenomes(nil, 1)
		g.Nts = make([][]byte, 1)
		g.Nts[0] = nts

		g.Names[0] = k
		fname := filepath.Join(ROOT, "All", fmt.Sprintf("%s.fasta", k))
		g.Save(k, fname, 0)
		fmt.Println("Wrote", fname)
	}
}

func main() {
	accessions := LoadAll()
	bySpecies := BySpecies(accessions)
	WhatsMissing(bySpecies)
	Assemble(bySpecies)
}
