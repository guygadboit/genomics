package main

import (
	"fmt"
	"genomics/genomes"
)

type Allele struct {
	Length int
	Pos    int
	Nts    map[string][]int
}

func FindAlleles(g *genomes.Genomes, length int) []Allele {
	ret := make([]Allele, 0)
	for i := 0; i < g.Length()-length; i++ {
		a := Allele{length, i, make(map[string][]int)}
		for j := 0; j < g.NumGenomes(); j++ {
			nts := string(g.Nts[j][i : i+length])
			_, there := a.Nts[nts]
			if !there {
				a.Nts[nts] = make([]int, 0)
			}
			a.Nts[nts] = append(a.Nts[nts], j)
		}
		ret = append(ret, a)
	}
	return ret
}

func (a *Allele) IsUnique(which int) (bool, string) {
	for k, v := range a.Nts {
		if len(v) == 1 && v[0] == which {
			return true, k
		}
	}
	return false, ""
}

func (a *Allele) Majority() string {
	ret := ""
	best := -1
	for k, v := range a.Nts {
		if len(v) > best {
			best = v
			ret = k
		}
	}
	return ret
}

// In all the places where which has a unique allele, show what the others have
func ShowUnique(alleles []Allele, which int) {
	for _, a := range alleles {
		if unique, us := a.IsUnique(which); unique {
			for k, v := range a.Nts {
				if k == us {
					continue
				}
				fmt.Printf("%d:%s %s: %d\n", a.Pos, us, k, len(v))
			}
			fmt.Println()
		}
	}
}

func MakeAncestor(g *genomes.Genomes, alleles []Allele, which int) {
	for _, a := range alleles {
		if a.IsUnique(which) {
			majority := Majority(a)
		}
	}
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	alleles := FindAlleles(g, 1)
	ShowUnique(alleles, 0)
}
