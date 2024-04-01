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
	return fmt.Sprintf("%s:%c%d%c", name, m.A, oPos/3+1, m.B)
}

func (m NtMut) ToString() string {
	var silence string
	if m.Silent {
		silence = "*"
	} else {
		silence = ""
	}
	return fmt.Sprintf("%c%d%c%s", m.A, m.Pos+1, m.B, silence)
}

type Comparison struct {
	Muts       []Mut
	NtMuts     []NtMut
	Insertions []int
	Deletions  []int
	Genomes    *genomes.Genomes
	A, B       int
}

// Returns silent and non-silent
func (c *Comparison) SilentCount() (S int, NS int) {
	for _, mut := range c.NtMuts {
		if mut.Silent {
			S++
		} else {
			NS++
		}
	}
	return
}

func (c *Comparison) Init(g *genomes.Genomes, a, b int) {
	c.Muts = make([]Mut, 0)
	c.NtMuts = make([]NtMut, 0)
	c.Insertions = make([]int, 0)
	c.Deletions = make([]int, 0)
	c.Genomes = g
	c.A, c.B = a, b
}

func (c *Comparison) Summary() {
	g := c.Genomes
	fmt.Printf("%d (%s) vs %d (%s)\n", c.A, g.Names[c.A], c.B, g.Names[c.B])
	fmt.Println("Amino acid changes")

	for _, mut := range c.Muts {
		fmt.Println(mut.ToString(c.Genomes.Orfs))
	}
	fmt.Printf("%d amino acids changed\n", len(c.Muts))

	fmt.Println("\nNucleotide changes")
	for _, mut := range c.NtMuts {
		fmt.Println(mut.ToString())
	}

	S, N := c.SilentCount()
	fmt.Printf("Silent: %d Non-Silent: %d Total: %d\n", S, N, S+N)
	fmt.Printf("Insertions: %d Deletions: %d\n",
		len(c.Insertions), len(c.Deletions))
}

func (c *Comparison) OneLineSummary() {
	g := c.Genomes
	fmt.Printf("%d-%d %s: ", c.A, c.B, g.Names[c.B])

	for _, mut := range c.Muts {
		fmt.Printf("%s ", mut.ToString(c.Genomes.Orfs))
	}

	for _, mut := range c.NtMuts {
		fmt.Printf("%s ", mut.ToString())
	}
	S, N := c.SilentCount()
	fmt.Printf("Silent: %d Non-Silent: %d Total: %d", S, N, S+N)
}

func compare(g *genomes.Genomes, a, b int) Comparison {
	var ret Comparison
	ret.Init(g, a, b)

	for _, aCodon := range genomes.Translate(g, a) {
		pos := aCodon.Pos
		if pos+3 >= len(g.Nts[b]) {
			break
		}
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
	var oneLine bool

	flag.StringVar(&orfName, "orfs", "", "ORFs")
	flag.StringVar(&include, "i", "", "Genomes to include (unset means all)")
	flag.StringVar(&saveTrans, "s", "", "Save Translation to file")
	flag.BoolVar(&oneLine, "1", false, "One line output")
	flag.Parse()

	var g *genomes.Genomes
	for i, fname := range flag.Args() {
		if i == 0 {
			g = genomes.LoadGenomes(fname, orfName, false)
			g.RemoveGaps()
		} else {
			g2 := genomes.LoadGenomes(fname, orfName, false)
			err := g.AlignCombine(g2)
			if err != nil {
				log.Print(err)
			}
		}
	}

	var which []int
	if include != "" {
		which = utils.ParseInts(include, ",")
	} else {
		which = make([]int, g.NumGenomes())
		for i := 0; i < g.NumGenomes(); i++ {
			which[i] = i
		}
	}

	for _, w := range which[1:] {
		c := compare(g, which[0], w)
		if oneLine {
			c.OneLineSummary()
		} else {
			c.Summary()
		}
		fmt.Println()
	}

	if len(which) == 2 && saveTrans != "" {
		g.SaveWithTranslation(saveTrans, nil, which[0], which[1])
		fmt.Printf("Wrote %s\n", saveTrans)
	}
}
