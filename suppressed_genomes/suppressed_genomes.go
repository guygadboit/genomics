package main

import (
	"fmt"
	"genomics/genomes"
	"log"
	"os"
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

func (a *Accession) MakeFname() string {
	ret := a.species
	if a.gene != Combined {
		ret = fmt.Sprintf("%s-%s.fasta", ret, GeneToString(a.gene))
	} else {
		ret += ".fasta"
	}
	return ret
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

func MakeIndex(accessions []Accession,
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

/*
Separate them out into Spikes, ORF8s etc and put each in its own directory
*/
func SaveByGene(byGene map[string][]Accession) {
	for k, v := range byGene {
		dir := filepath.Join(ROOT, k)
		err := os.MkdirAll(dir, 0750)
		if err != nil {
			log.Fatal(err)
		}

		for _, acc := range v {
			fname := filepath.Join(dir, acc.MakeFname())
			acc.genome.Save(acc.genome.Names[0], fname, 0)
			fmt.Printf("Wrote %s\n", fname)
		}
	}
}

/*
Merge ORF8, RdRp and Spike and save the most complete genome we have for each
species
*/
func SaveBySpecies(bySpecies map[string][]Accession) {
	for k, v := range bySpecies {
		g := genomes.NewGenomes(nil, 1)
		g.Nts[0] = make([]byte, 0)

		for _, acc := range v {
			g.Nts[0] = append(g.Nts[0], acc.genome.Nts[0]...)
		}

		fname := filepath.Join(ROOT, "All", fmt.Sprintf("%s.fasta", k))
		g.Save(k, fname, 0)
		fmt.Println("Wrote", fname)
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

	bySpecies := MakeIndex(accessions, func(a *Accession) string {
		return a.species
	})
	SaveBySpecies(bySpecies)

}
