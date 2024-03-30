package main

import (
	"fmt"
	"log"
	"flag"
	"genomics/utils"
	"genomics/genomes"
)

// Represents an AA change
type Mut struct {
	A, B	byte
	Pos		int
}

// Represents a NT change
type NtMut struct {
	Mut
	Silent	bool
}

func (m Mut) ToString(orfs genomes.Orfs) string {
	orfI, oPos, err := orfs.GetOrfRelative(m.Pos)
	if err != nil {
		return fmt.Sprintf("Not in ORF")
	}
	name := orfs[orfI].Name
	return fmt.Sprintf("%s:%c%d%c", name, m.A, oPos/3+1, m.B)
}

func (m NtMut) ToString() string {
	var silence byte
	if m.Silent {
		silence = 'S'
	} else {
		silence = 'N'
	}
	return fmt.Sprintf("%c%d%c %c", m.A, m.Pos+1, m.B, silence)
}

type Comparison struct {
	Muts	[]Mut
	NtMuts	[]NtMut
	Orfs	genomes.Orfs
}

func (c *Comparison) Init(orfs genomes.Orfs) {
	c.Muts = make([]Mut, 0)
	c.NtMuts = make([]NtMut, 0)
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
}

func compare(g *genomes.Genomes, a, b int) Comparison {
	var ret Comparison
	ret.Init(g.Orfs)

	for _, aCodon := range genomes.Translate(g, a) {
		pos := aCodon.Pos
		bNts := string(g.Nts[b][pos:pos+3])

		if aCodon.Nts == bNts {
			continue
		}

		var bCodon genomes.Codon
		bCodon.Init(pos, bNts)
		silent := aCodon.Aa == bCodon.Aa

		for j, aNt := range []byte(aCodon.Nts) {
			bNt := bCodon.Nts[j]
			if aNt != bNt {
				ret.NtMuts = append(ret.NtMuts,
					NtMut{Mut{aNt, bNt, pos+j}, silent})
			}
		}

		if !silent {
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
