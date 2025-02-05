package comparison

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
	"slices"
)

// Represents an AA change
type Mut struct {
	A, B byte
	Pos  int
}

// Represents a nt change
type NtMut struct {
	Mut
	Silence utils.Silence
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
	case utils.SILENT:
		silence = "*"
	case utils.NON_SILENT:
		silence = ""
	case utils.NOT_IN_ORF:
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
		case utils.SILENT:
			S++
		case utils.NON_SILENT:
			NS++
		case utils.NOT_IN_ORF:
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
		if mut.Silence == utils.SILENT {
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

func (c *Comparison) Summary(showIndels bool) {
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

	fmt.Printf("\n%d nucleotides changed\n", len(c.NtMuts))
	for _, mut := range c.NtMuts {
		fmt.Println(mut.ToString())
	}

	if showIndels {
		fmt.Println("Insertions (1-based)")
		for _, ins := range c.Insertions {
			fmt.Println(ins + 1)
		}

		fmt.Println("Deletions (1-based)")
		for _, del := range c.Deletions {
			fmt.Println(del + 1)
		}
	}

	S, N, NO := c.SilentCount()
    ratio := float64(N)/float64(S)
    fmt.Printf("Silent: %d Non-Silent: %d Non-Orf: %d N/S: %.2f Total: %d\n",
		S, N, NO, ratio, S+N+NO)
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

type Transition struct {
	from, to byte
}

type TransitionCount struct {
	t     Transition
	count int
}

func (c *Comparison) ShowTransitions() {
	makeTransitions := func(silentOnly bool) map[Transition]int {
		transitions := make(map[Transition]int)
		for _, ntm := range c.NtMuts {
			if !silentOnly || ntm.Silence != utils.NON_SILENT {
				t := Transition{ntm.A, ntm.B}
				transitions[t]++
			}
		}
		return transitions
	}

	makeCounts := func(t map[Transition]int) ([]TransitionCount, int) {
		ret := make([]TransitionCount, 0, len(t))
		total := 0
		for k, v := range t {
			ret = append(ret, TransitionCount{k, v})
			total += v
		}

		slices.SortFunc(ret, func(a, b TransitionCount) int {
			if a.count < b.count {
				return 1
			}
			if a.count > b.count {
				return -1
			}
			return 0
		})
		return ret, total
	}

	showCounts := func(counts []TransitionCount, total int) {
		for _, c := range counts {
			rate := float64(c.count) / float64(total)
			fmt.Printf("%c->%c %d/%d %.2f\n", c.t.from,
				c.t.to, c.count, total, rate)
		}
	}

	fmt.Println("Transitions (all)")
	counts, total := makeCounts(makeTransitions(false))
	showCounts(counts, total)

	fmt.Println("Transitions (silent only)")
	counts, total = makeCounts(makeTransitions(true))
	showCounts(counts, total)
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
		if m.Silence == utils.NON_SILENT {
			nonSilent[m.Pos] = true
		} else {
			silent[m.Pos] = true
		}
	}

	insertions := utils.ToSet(c.Insertions)
	deletions := utils.ToSet(c.Deletions)

	g := c.Genomes

	var s, ns, ins, del int
	var nsRatio float64

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
		if s > 0 {
			nsRatio = float64(ns) / float64(s)
		}

		// The ratio goes in the data file but isn't plotted by default
		fmt.Fprintf(w, "%d %d %d %d %.4f\n", s, ns, ins, del, nsRatio)
	}

	w.Flush()
}

func Compare(g *genomes.Genomes, a, b int) Comparison {
	var ret Comparison
	ret.Init(g, a, b, false)

	handleNtMut := func(aNt, bNt byte, silence utils.Silence, pos int) bool {
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
		var silence utils.Silence
		var isMut bool

		isMut = aCodon.Aa != '-' && bCodon.Aa != '-'
		if isMut {
			if aCodon.Aa == bCodon.Aa {
				silence = utils.SILENT
			} else {
				silence = utils.NON_SILENT
			}
		}

		for j, aNt := range []byte(aCodon.Nts) {
			bNt := bCodon.Nts[j]

			if !handleNtMut(aNt, bNt, silence, pos+j) {
				continue
			}
		}

		if isMut && silence == utils.NON_SILENT {
			ret.Muts = append(ret.Muts, Mut{aCodon.Aa, bCodon.Aa, pos})
		}
	}

	// And now any outside them.
	start := 0
	for _, orf := range g.Orfs {
		for i := start; i < orf.Start; i++ {
			aNt, bNt := g.Nts[a][i], g.Nts[b][i]
			handleNtMut(aNt, bNt, utils.NOT_IN_ORF, i)

		}
		start = orf.End
	}

	slices.Sort(ret.Insertions)
	slices.Sort(ret.Deletions)

	return ret
}

// Use this when the genomes are already translated
func CompareProtein(g *genomes.Genomes, a, b int) Comparison {
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
