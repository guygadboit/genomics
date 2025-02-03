package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
)

type Allele struct {
	Pos int
	Nts map[string][]int
}

type Alleles struct {
	Length     int
	NumGenomes int
	Data       []Allele
}

func NewAlleles(g *genomes.Genomes, length int) *Alleles {
	return &Alleles{length, g.NumGenomes(), make([]Allele, 0)}
}

func FindAlleles(g *genomes.Genomes, length int) *Alleles {
	ret := NewAlleles(g, length)
	for i := 0; i < g.Length()+1-length; i++ {
		a := Allele{i, make(map[string][]int)}
		for j := 0; j < g.NumGenomes(); j++ {
			nts := string(g.Nts[j][i : i+length])

			if length == 1 && !utils.IsRegularNt(nts[0]) {
				continue
			}

			_, there := a.Nts[nts]
			if !there {
				a.Nts[nts] = make([]int, 0)
			}
			a.Nts[nts] = append(a.Nts[nts], j)
		}
		ret.Data = append(ret.Data, a)
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
			best = len(v)
			ret = k
		}
	}
	return ret
}

// In all the places where which has a unique allele, show what the others have
func (al *Alleles) ShowUnique(which int) {
	fmt.Printf("Unique alleles\n")
	for _, a := range al.Data {
		if unique, us := a.IsUnique(which); unique {
			got := false
			for k, v := range a.Nts {
				if k == us {
					continue
				}
				prop := float64(len(v)) / float64(al.NumGenomes)
				fmt.Printf("%d%s %d/%d have %s (%.2f)\n",
					a.Pos+1, us, len(v), al.NumGenomes, k, prop)
				got = true
			}
			if got {
				fmt.Println()
			}
		}
	}
}

func (al *Alleles) MakeAncestor(g *genomes.Genomes,
	which int) *genomes.Genomes {
	if al.Length != 1 {
		log.Fatal("This only works if your window is 1")
	}

	changes := 0
	ret := genomes.NewGenomes(g.Orfs, 1)
	ret.Nts[0] = make([]byte, len(al.Data))

	fmt.Printf("Mutations being applied:\n")
	for i, a := range al.Data {
		nt := g.Nts[which][i]
		if unique, _ := a.IsUnique(which); unique {
			majority := a.Majority()[0]
			if nt == majority {
				fmt.Printf("Already have %d%c\n", i+1, nt)
			} else {
				changes++
				fmt.Printf("%c%d%c\n", nt, i+1, majority)
			}
			nt = majority
		}
		ret.Nts[0][i] = nt
	}
	fmt.Printf("%d changes\n", changes)
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	alleles := FindAlleles(g, 1)
	alleles.ShowUnique(0)
	ancestor := alleles.MakeAncestor(g, 0)
	ancestor.Save("Ancestor", "ancestor.fasta", 0)
	fmt.Printf("Wrote ancestor.fasta\n")
}
