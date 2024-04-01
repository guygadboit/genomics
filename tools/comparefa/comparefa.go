package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
)

// Represents an AA change
type Mut struct {
	A, B byte
	Pos  int
}

// Represents a nt change
type NtMut struct {
	Mut
	Silent bool
}

func (m Mut) ToString(orfs genomes.Orfs) string {
	orfI, oPos, err := orfs.GetOrfRelative(m.Pos)
	if err != nil {
		return fmt.Sprintf("Not in ORF")
	}
	name := orfs[orfI].Name
	return fmt.Sprintf("%s\t\t%c%d%c", name, m.A, oPos/3+1, m.B)
}

func (m NtMut) ToString() string {
	var silence string
	if m.Silent {
		silence = "S"
	} else {
		silence = "NS"
	}
	return fmt.Sprintf("%c%d%c %s", m.A, m.Pos+1, m.B, silence)
}

type Comparison struct {
	Muts     []Mut
	NtMuts   []NtMut
	Insertions	[]int
	Deletions	[]int
	Orfs     genomes.Orfs
}

func (c *Comparison) Init(orfs genomes.Orfs) {
	c.Muts = make([]Mut, 0)
	c.NtMuts = make([]NtMut, 0)
	c.Insertions = make([]int, 0)
	c.Deletions = make([]int, 0)
	c.Orfs = orfs
}

func (c Comparison) Summary() {
	fmt.Println("Amino acid Changes")

	for _, mut := range c.Muts {
		fmt.Println(mut.ToString(c.Orfs))
	}
	fmt.Printf("%d amino acids changed\n", len(c.Muts))

	var S, N int
	fmt.Println("\nNucleotide Changes")
	for _, mut := range c.NtMuts {
		fmt.Println(mut.ToString())
		if mut.Silent {
			S++
		} else {
			N++
		}
	}
	fmt.Printf("Silent: %d Non-Silent: %d Total: %d\n", S, N, S+N)
	fmt.Printf("Insertions: %d Deletions: %d\n",
		len(c.Insertions), len(c.Deletions))
}

func compare(g *genomes.Genomes, a, b int) Comparison {
	var ret Comparison
	ret.Init(g.Orfs)

	for _, aCodon := range genomes.Translate(g, a) {
		pos := aCodon.Pos
		bNts := string(g.Nts[b][pos : pos+3])

		if aCodon.Nts == bNts {
			continue
		}

		var bCodon genomes.Codon
		bCodon.Init(pos, bNts)
		var silent bool
		var isMut bool

		isMut = aCodon.Aa != '-' && bCodon.Aa != '-'
		silent = isMut && aCodon.Aa == bCodon.Aa

		for j, aNt := range []byte(aCodon.Nts) {
			bNt := bCodon.Nts[j]
			if aNt == bNt {
				continue
			}

			if aNt == 'N' || bNt == 'N' {
				continue
			}

			if aNt == '-' {
				ret.Insertions = append(ret.Insertions, pos)
			} else if bNt == '-' {
				ret.Deletions = append(ret.Deletions, pos)
			} else {
				ret.NtMuts = append(ret.NtMuts,
					NtMut{Mut{aNt, bNt, pos + j}, silent})
			}
		}

		if isMut && !silent {
			ret.Muts = append(ret.Muts, Mut{aCodon.Aa, bCodon.Aa, pos})
		}
	}
	return ret
}

func main() {
	var orfName string
	var include string
	var saveTrans string

	flag.StringVar(&orfName, "orfs", "", "ORFs")
	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.StringVar(&saveTrans, "save", "", "Save Translation to file")
	flag.Parse()

	var which []int
	if include != "" {
		which = utils.ParseInts(include, ",")
		if len(which) != 2 {
			log.Fatal("Only two genomes for now")
		}
	} else {
		which = []int{0, 1}
	}

	g := genomes.LoadGenomes(flag.Arg(0), orfName, false)
	c := compare(g, which[0], which[1])
	c.Summary()

	if saveTrans != "" {
		g.SaveWithTranslation(saveTrans, nil, which[0], which[1])
		fmt.Printf("Wrote %s\n", saveTrans)
	}
}
