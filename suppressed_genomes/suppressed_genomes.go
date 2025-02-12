package main

import (
	"fmt"
	"genomics/genomes"
	"log"
	"path/filepath"
	"regexp"
	"slices"
	"strings"
)

var ROOT string = "/fs/f/genomes/viruses/suppressed_genomes/"

type Gene int

const (
	RdRp Gene = iota
	S
	ORF8
	Combined // Means all of them together
)

var GENES []string = []string{"RdRP", "S", "ORF8", "Combined"}

func GeneFromString(g string) Gene {
	index := slices.Index(GENES, g)
	if index != -1 {
		return Gene(index)
	}
	if strings.Contains(g, "RdRp") {
		return RdRp
	}
	if strings.Contains(g, "spike") {
		return S
	}
	log.Fatalf("Unrecognized gene %s\n", g)
	return -1
}

func GeneToString(g Gene) string {
	genes := []string{"RdRp", "S", "ORF8"}
	return genes[g]
}

type Accession struct {
	fname  string
	genome *genomes.Genomes
	gene   Gene

	// The genome name minus the bit that tells you it's ORF8 or whatever. So
	// we can match them up.
	species  string
	location string
}

func (a *Accession) ToString() string {
	return fmt.Sprintf("%s %s %s", a.species, GeneToString(a.gene), a.location)
}

func LoadAll() []Accession {
	pat := regexp.MustCompile(`.*coronavirus strain ([^\s]+) (.+) gene.*$`)
	locPat := regexp.MustCompile(`^.*_(.*)`)

	ret := make([]Accession, 0)
	fullPath := filepath.Join(ROOT, "Downloads", "MH*.fasta.gz")
	matches, _ := filepath.Glob(fullPath)
	for _, m := range matches {
		g := genomes.LoadGenomes(m, "", false)

		name := g.Names[0]
		fmt.Println(name)
		matches := pat.FindSubmatch([]byte(name))
		species := string(matches[1])
		gene := string(matches[2])

		matches = locPat.FindSubmatch(matches[1])
		location := string(matches[1])

		acc := Accession{m, g, GeneFromString(gene), species, location}
		ret = append(ret, acc)
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

func MakeIndex(accessions[] Accession,
	keyFn func(a *Accession) string) map[string][]Accession {
	ret := make(map[string][]Accession)
	for _, acc := range accessions {
		key := keyFn(&acc)
		_, there := ret[key]
		if !there {
			ret[key] = make([]Accession, 0)
		}
		ret[key] = append(ret[key], acc)
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

// FIXME this is no good. We want one BySpecies, which puts the ORFs together,
// and then another ByLocation, which is a completeish genome for each species
// in a fasta file (that could in theory then be aligned)
func Assemble(byMap map[string][]Accession, merge bool) {
	for k, v := range byMap {
		numGenomes := 1
		if !merge {
			numGenomes = len(v)
		}

		g := genomes.NewGenomes(nil, numGenomes)
		g.Nts = make([][]byte, numGenomes)

		if merge {
			g.Nts[0] = make([]byte, 0)
		}

		for i, acc := range v {
			if merge {
				g.Nts[0] = append(g.Nts[0], acc.genome.Nts[0]...)
			} else {
				g.Nts[i] = acc.genome.Nts[0]
				g.Names[i] = acc.genome.Names[0]
			}
		}

		fname := filepath.Join(ROOT, "All", fmt.Sprintf("%s.fasta", k))
		g.SaveMulti(fname)
		fmt.Println("Wrote", fname)
	}
}

func SaveByGene(byGene map[string][]Accession) {
	for k, v := range byGene {
		dir := filepath.Join(ROOT, k)
		err := os.MkdirAll(dir, 0750)
		if err != nil {
			log.Fatal(err)
		}

		// FIXME you are here

	}
}

func main() {
	accessions := LoadAll()
	for _, acc := range accessions {
		fmt.Println(acc.ToString())
	}

	byGene := MakeIndex(accessions, func(a *Accession) string {
		return GeneToString(a.gene)
	})

	SaveByGene(byGene)

}
