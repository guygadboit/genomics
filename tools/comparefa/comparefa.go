package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
)

// Represents an AA change
type Mut struct {
	A, B byte
	Pos  int
}

type Silence int

const (
	SILENT Silence = iota
	NON_SILENT
	NOT_IN_ORF
)

// Represents a nt change
type NtMut struct {
	Mut
	Silence Silence
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

	switch m.Silence {
	case SILENT:
		silence = "*"
	case NON_SILENT:
		silence = ""
	case NOT_IN_ORF:
		silence = "@"
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
	IsProtein  bool
}

// Returns silent, non-silent and non-orf
func (c *Comparison) SilentCount() (S int, NS int, NO int) {
	for _, mut := range c.NtMuts {
		switch mut.Silence {
		case SILENT:
			S++
		case NON_SILENT:
			NS++
		case NOT_IN_ORF:
			NO++
		}
	}
	return
}

func (c *Comparison) Init(g *genomes.Genomes, a, b int, protein bool) {
	c.Muts = make([]Mut, 0)
	c.NtMuts = make([]NtMut, 0)
	c.Insertions = make([]int, 0)
	c.Deletions = make([]int, 0)
	c.Genomes = g
	c.A, c.B = a, b
	c.IsProtein = protein
}

func (c *Comparison) SilentSummary() {
	for _, mut := range c.NtMuts {
		if mut.Silence == SILENT {
			fmt.Println(mut.ToString())
		}
	}
}

func (c *Comparison) proteinSummary() {
	g := c.Genomes
	fmt.Printf("%d (%s) vs %d (%s)\n", c.A, g.Names[c.A], c.B, g.Names[c.B])
	fmt.Println("Amino acid changes")

	for _, mut := range c.Muts {
		fmt.Printf("%c%d%c\n", mut.A, mut.Pos, mut.B)
	}
	fmt.Printf("%d amino acids changed\n", len(c.Muts))

	numMuts := float64(len(c.Muts))
	total := float64(g.Length())
	fmt.Printf("AA similarity: %.2f%%\n", ((total-numMuts)*100)/total)
}

func (c *Comparison) Summary() {
	if c.IsProtein {
		c.proteinSummary()
		return
	}
	g := c.Genomes
	fmt.Printf("%d (%s) vs %d (%s)\n", c.A, g.Names[c.A], c.B, g.Names[c.B])
	fmt.Println("Amino acid changes")

	for _, mut := range c.Muts {
		fmt.Printf("%s (%d)\n", mut.ToString(c.Genomes.Orfs), mut.Pos+1)
	}
	fmt.Printf("%d amino acids changed\n", len(c.Muts))

	fmt.Println("\nNucleotide changes")
	for _, mut := range c.NtMuts {
		fmt.Println(mut.ToString())
	}

	S, N, NO := c.SilentCount()
	fmt.Printf("Silent: %d Non-Silent: %d Non-Orf: %d Total: %d\n",
		S, N, NO, S+N+NO)
	fmt.Printf("Insertions: %d Deletions: %d\n",
		len(c.Insertions), len(c.Deletions))

	numMuts := float64(len(c.NtMuts))
	total := float64(g.Length())
	fmt.Printf("Nucleotide similarity: %.2f%%\n", ((total-numMuts)*100)/total)
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
	S, N, NO := c.SilentCount()
	fmt.Printf("Silent: %d Non-Silent: %d Non-Orf: %d Total: %d\n",
		S, N, NO, S+N+NO)
}

type PosSet map[int]bool

func (c *Comparison) GraphData(fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)

	silent := make(PosSet)
	nonSilent := make(PosSet)
	for _, m := range c.NtMuts {
		if m.Silence == NON_SILENT {
			nonSilent[m.Pos] = true
		} else {
			silent[m.Pos] = true
		}
	}

	insertions := utils.ToSet(c.Insertions)
	deletions := utils.ToSet(c.Deletions)

	g := c.Genomes

	var s, ns, ins, del int
	for i := 0; i < g.Length(); i++ {
		if silent[i] {
			s++
		} else if nonSilent[i] {
			ns++
		} else if insertions[i] {
			ins++
		} else if deletions[i] {
			del++
		}
		fmt.Fprintf(w, "%d %d %d %d\n", s, ns, ins, del)
	}

	w.Flush()
}

func compare(g *genomes.Genomes, a, b int) Comparison {
	var ret Comparison
	ret.Init(g, a, b, false)

	handleNtMut := func(aNt, bNt byte, silence Silence, pos int) bool {
		if aNt == bNt {
			return false
		}
		if aNt == 'N' || bNt == 'N' {
			return false
		}
		if aNt == '-' {
			ret.Insertions = append(ret.Insertions, pos)
		} else if bNt == '-' {
			ret.Deletions = append(ret.Deletions, pos)
		} else {
			ret.NtMuts = append(ret.NtMuts,
				NtMut{Mut{aNt, bNt, pos}, silence})
		}
		return true
	}

	// First find the mutations in ORFs
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
		var silence Silence
		var isMut bool

		isMut = aCodon.Aa != '-' && bCodon.Aa != '-'
		if isMut {
			if aCodon.Aa == bCodon.Aa {
				silence = SILENT
			} else {
				silence = NON_SILENT
			}
		}

		for j, aNt := range []byte(aCodon.Nts) {
			bNt := bCodon.Nts[j]

			if !handleNtMut(aNt, bNt, silence, pos+j) {
				continue
			}
		}

		if isMut && silence == NON_SILENT {
			ret.Muts = append(ret.Muts, Mut{aCodon.Aa, bCodon.Aa, pos})
		}
	}

	// And now any outside them.
	start := 0
	for _, orf := range g.Orfs {
		for i := start; i < orf.Start; i++ {
			aNt, bNt := g.Nts[a][i], g.Nts[b][i]
			handleNtMut(aNt, bNt, NOT_IN_ORF, i)

		}
		start = orf.End
	}

	return ret
}

// Use this when the genomes are already translated
func compareProtein(g *genomes.Genomes, a, b int) Comparison {
	var ret Comparison
	ret.Init(g, a, b, true)

	aAas, bAas := g.Nts[a], g.Nts[b]

	for i := 0; i < g.Length(); i++ {
		aAa, bAa := aAas[i], bAas[i]
		if aAa == bAa {
			continue
		}
		if aAa == '-' {
			ret.Insertions = append(ret.Insertions, i)
		} else if bAa == '-' {
			ret.Deletions = append(ret.Deletions, i)
		} else {
			ret.Muts = append(ret.Muts, Mut{aAa, bAa, i})
		}
	}
	return ret
}

func main() {
	var orfName string
	var include string
	var saveTrans string
	var oneLine bool
	var keepGaps bool
	var graphData bool
	var silentOnly bool
	var protein bool
	var highlightString string

	flag.StringVar(&orfName, "orfs", "", "ORFs")
	flag.StringVar(&include, "i", "", "Genomes to include (unset means all)")
	flag.StringVar(&saveTrans, "s", "", "Save translation to file")
	flag.BoolVar(&oneLine, "1", false, "One line output")
	flag.BoolVar(&keepGaps, "gaps", false, "Keep gaps in first genome")
	flag.BoolVar(&graphData, "g", false, "Graph data")
	flag.BoolVar(&silentOnly, "silent", false, "Silent only")
	flag.BoolVar(&protein, "p", false, "input is already translated")
	flag.StringVar(&highlightString, "highlights",
		"", "1-based positions to highlight separated with ,")
	flag.Parse()

	var g *genomes.Genomes
	for i, fname := range flag.Args() {
		if i == 0 {
			g = genomes.LoadGenomes(fname, orfName, false)
			if !keepGaps {
				g.RemoveGaps()
			}
		} else {
			g2 := genomes.LoadGenomes(fname, orfName, false)
			err := g.AlignCombine(g2)
			if err != nil {
				log.Print(err)
			}
		}
	}

	err := g.CheckOrfs()
	if err != nil {
		fmt.Fprintln(os.Stderr, "Maybe need -gaps?")
		log.Fatal(err)
	}

	var which []int
	if include != "" {
		which = utils.ParseInts(include, ",")
		for _, v := range which {
			if v < 0 || v > g.NumGenomes()-1 {
				log.Fatalf("Invalid index %d", v)
			}
		}
	} else {
		which = make([]int, g.NumGenomes())
		for i := 0; i < g.NumGenomes(); i++ {
			which[i] = i
		}
	}

	for _, w := range which[1:] {
		var c Comparison
		if protein {
			c = compareProtein(g, which[0], w)
		} else {
			c = compare(g, which[0], w)
		}
		if oneLine {
			c.OneLineSummary()
		} else if silentOnly {
			c.SilentSummary()
		} else {
			c.Summary()
		}
		if len(which) == 2 && graphData {
			fname := "cumulative-muts.txt"
			c.GraphData(fname)
			fmt.Printf("Wrote %s\n", fname)
		}
	}

	if len(which) == 2 && saveTrans != "" {
		highlights := genomes.ParseHighlights(highlightString, ",", true, 'v')
		g.SaveWithTranslation(saveTrans, highlights, which[0], which[1])
		fmt.Printf("Wrote %s\n", saveTrans)
	}
}
