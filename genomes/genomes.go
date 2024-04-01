package genomes

import (
	"bufio"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"slices"
	"strings"
	"errors"
)

/*
Represents a collection of aligned genomes (usually two) with one genome
per row. The orfs "belong" to the first one in the set. We also use this
for a single genome.
*/
type Genomes struct {
	Nts   [][]byte
	Names []string
	Orfs  Orfs
}

func NewGenomes(orfs Orfs, numGenomes int) *Genomes {
	return &Genomes{make([][]byte, numGenomes),
		make([]string, numGenomes), orfs}
}

/*
Load genomes, which might be a fasta file containing a single genome, or
one containing a few of them in an alignment. Be a bit careful when working
with alignments since there may be '-' in there. You might want to call
RemoveGaps before going any futher. If merge load everything into one
genome even if there are multiple > lines (mammal genomes tend to be like
that). orfsName can be empty string if you don't have any ORFs which is
fine if you don't plan on doing any translation.
*/
func LoadGenomes(fname string, orfsName string, merge bool) *Genomes {
	var orfs Orfs

	if orfsName != "" {
		orfs = LoadOrfs(orfsName)
	}

	ret := NewGenomes(orfs, 0)

	fp := utils.NewFileReader(fname)
	defer fp.Close()

	currentRow := make([]byte, 0)
loop:
	// If working with huge genomes uncomment this for faster debugging!
	// for i := 0; i < 1000; i++ {
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)

		if strings.HasPrefix(line, ">") {
			/*
				fields := strings.Fields(line[1:])
				name := fields[0]
			*/
			name := line[1:]
			ret.Names = append(ret.Names, name)
			if !merge && len(currentRow) > 0 {
				ret.Nts = append(ret.Nts, currentRow)
				currentRow = make([]byte, 0)
			}
			continue
		}
		line = strings.ToUpper(line)
		currentRow = append(currentRow, []byte(line)...)
	}
	ret.Nts = append(ret.Nts, currentRow)

	if len(ret.Orfs) == 0 {
		ret.Orfs = []Orf{{0, ret.Length(), ""}}
	}

	return ret
}

func (g *Genomes) PrintSummary() {
	for i, n := range g.Names {
		fmt.Printf("%d: %s\n", i, n)
	}
}

func (g *Genomes) Save(name, fname string, which int) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	fmt.Fprintf(fp, ">%s\n", name)

	nts := g.Nts[which]
	utils.Wrap(fp, nts)
	fp.Flush()
	return nil
}

/*
Make one fasta file with all the nt sequences in it and their names.
*/
func (g *Genomes) SaveMulti(fname string) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)

	for i := 0; i < g.NumGenomes(); i++ {
		fmt.Fprintf(fp, ">%s\n", g.Names[i])
		utils.Wrap(fp, g.Nts[i])
	}
	fp.Flush()
	return nil
}

func (g *Genomes) SaveSelected(fname string, which ...int) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)

	for _, i := range which {
		fmt.Fprintf(fp, ">%s\n", g.Names[i])
		nts := g.Nts[i]
		utils.Wrap(fp, nts)
	}

	fp.Flush()
	return nil
}

func (g *Genomes) Clone() *Genomes {
	ret := NewGenomes(g.Orfs, g.NumGenomes())
	for i := 0; i < len(g.Nts); i++ {
		ret.Nts[i] = make([]byte, len(g.Nts[i]))
		copy(ret.Nts[i], g.Nts[i])
	}
	return ret
}

// Make a shallow copy that matches the specified genomes
func (g *Genomes) Filter(which ...int) *Genomes {
	ret := NewGenomes(g.Orfs, len(which))

	for i, w := range which {
		ret.Nts[i] = g.Nts[w]
		ret.Names[i] = g.Names[w]
	}

	return ret
}

// Make the specified genomes into deep copies of themselves (because you want
// to mutate them for example)
func (g *Genomes) DeepCopy(which ...int) {
	n := g.Length()
	for _, w := range which {
		newNts := make([]byte, n)
		copy(newNts, g.Nts[w])
		g.Nts[w] = newNts
	}
}

/*
Assuming align is aligned with g, add it to g's own nts array, just doing a
shallow copy
*/
func (g *Genomes) Combine(other *Genomes) error {
	for i := 0; i < other.NumGenomes(); i++ {
		if len(other.Nts[i]) != len(g.Nts[0]) {
			return errors.New("Genomes aren't the same length")
		}
		g.Nts = append(g.Nts, other.Nts[i])
	}
	g.Names = append(g.Names, other.Names...)
	return nil
}

// Insert '-' into a where necessary to bring it up to the length of ref and
// maintaining a good alignment.
func handleDeletions(a []byte, ref []byte) {
	needed := len(ref) - len(a)

	// Would putting a gap this long in at pos be a good idea?
	tryGap := func(pos, length int) bool {
		overflow := pos+length+10 - len(a)
		if overflow > 0 {
			length -= overflow
		}

		var good int
		for i := pos; i < pos+10; i++ {
			if a[i] == ref[i+length] {
				good++
			}
		}

		return good >= 8
	}

	// Actually put a gap in
	insertGap := func(pos, length int) {
		// Weird Go idiom for inserting 
		a = append(a[:pos], append(make([]byte, length), a[pos:]...)...)
		for i := pos; i < pos+length; i++ {
			a[i] = '-'
		}
		fmt.Printf("Inserted %d at %d\n", length, pos)
	}

	for i := 0; i < len(ref); i++ {
		if i < len(a) {
			if a[i] == ref[i] {
				continue
			}
			for j := 1; j < needed; j++ {
				if tryGap(i, j) {
					insertGap(i, j)
					needed -= j
					i += j
					break
				}
			}
		}
	}

	for i := 0; i < needed; i++ {
		a = append(a, '-')
	}
}

/*
Quick and dirty alignment. We aren't sure yet what we actually need here so
we'll just start with padding.
*/
func (g *Genomes) AlignCombine(other *Genomes) error {
	for i := 0; i < other.NumGenomes(); i++ {
		n := len(other.Nts[i])

		// The other genome contains an insertion. For now just skip it, so we
		// can see how many there are.
		if n > len(g.Nts[0]) {
			log.Printf("%d contains %d insertions\n",
				i, len(other.Nts[i]) - len(g.Nts[0]))
			continue
		} else {
			handleDeletions(other.Nts[i], g.Nts[0])
		}

		g.Nts = append(g.Nts, other.Nts[i])
	}
	g.Names = append(g.Names, other.Names...)
	return nil
}

/*
These little functions make it a bit easier not to get confused about
which dimension is which.
*/
func (g *Genomes) NumGenomes() int {
	return len(g.Nts)
}

func (g *Genomes) Length() int {
	return len(g.Nts[0])
}

func (g *Genomes) Slice(which, start, end int) []byte {
	nts := g.Nts[which]

	if start < 0 {
		start = 0
	}
	if end > len(nts) {
		end = len(nts)
	}
	return nts[start:end]
}

/*
Whereever there is a gap in the first genome of an alignment, just remove that
column. You will need to do this for comparisons involving translations to work
as the Orfs don't take into account gaps, and neither do Environment
operations. Returns how many we dropped.
*/
func (g *Genomes) RemoveGaps() int {
	var ret int
	nts := g.Nts
	newNts := make([][]byte, g.NumGenomes())

	for i := 0; i < g.NumGenomes(); i++ {
		newNts[i] = make([]byte, 0, g.Length())
	}

	for i := 0; i < g.Length(); i++ {
		if nts[0][i] == '-' {
			ret++
			continue
		}

		for j := 0; j < g.NumGenomes(); j++ {
			newNts[j] = append(newNts[j], nts[j][i])
		}
	}

	g.Nts = newNts
	return ret
}

func (g *Genomes) HaveOrfs() bool {
	return len(g.Orfs) > 0
}

/*
What's the sequence similarity between the a'th and b'th genomes in an
alignment? Returns a number between 0 and 1 (so multiply by 100 if you want a
percentage)
*/
func (g *Genomes) SequenceSimilarity(a, b int) float64 {
	var same, total int
	for i := 0; i < g.Length(); i++ {
		// Only consider proper nts
		switch g.Nts[a][i] {
		case 'A':
			fallthrough
		case 'C':
			fallthrough
		case 'G':
			fallthrough
		case 'T':
			break
		default:
			continue
		}
		total++
		if g.Nts[a][i] == g.Nts[b][i] {
			same++
		}
	}
	return float64(same) / float64(total)
}

/*
What is the sequence similarity between (start, end]? If include is false, then
return the SS everywhere else instead
*/
func (g *Genomes) SubSequenceSimilarity(a, b int,
	start, end int, include bool) float64 {
	var same, total int
	for i := 0; i < g.Length(); i++ {
		// Only consider proper nts
		switch g.Nts[a][i] {
		case 'A':
			fallthrough
		case 'C':
			fallthrough
		case 'G':
			fallthrough
		case 'T':
			break
		default:
			continue
		}

		inRange := i >= start && i < end
		if inRange != include {
			continue
		}

		total++
		if g.Nts[a][i] == g.Nts[b][i] {
			same++
		}
	}
	return float64(same) / float64(total)
}

func shorten(s string, length int) string {
	words := strings.Split(s, " ")
	w := words[0]
	if len(w) > length {
		return w[0:length]
	} else {
		return w
	}
}

// Return a string with * where they're all the same or ' ' when not.
func (g *Genomes) compare(start, end int, which ...int) string {
	n := end - start
	ret := make([]byte, n)
outer:
	for i := 0; i < n; i++ {
		pos := i + start
		nts := make(map[byte]bool)
		for j, w := range which {
			if j > 0 {
				_, there := nts[g.Nts[w][pos]]
				if !there {
					ret[i] = ' '
					continue outer
				}
			}
			nts[g.Nts[w][pos]] = true
		}
		ret[i] = '*'
	}
	return string(ret)
}

// Mark a region with some character above the alignment.
type Highlight struct {
	Start int
	End   int
	Char  byte
}

func (g *Genomes) saveCluStyle(fname string,
	highlights []Highlight, withTranslation bool, which ...int) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()
	fp := bufio.NewWriter(fd)

	if which[0] == -1 {
		which = make([]int, g.NumGenomes())
		for i := 0; i < g.NumGenomes(); i++ {
			which[i] = i
		}
	}

	names := make([]string, len(which))
	trans := make([]Translation, len(which))
	var index, nextIndex int
	for i, w := range which {
		names[i] = shorten(g.Names[w], 12)
		if withTranslation {
			trans[i] = Translate(g, w)
		}
	}

	slices.SortFunc(highlights, func(a, b Highlight) int {
		return a.Start - b.Start
	})
	var highlightIt int

	for i := 0; i < g.Length(); i += 60 {
		n := g.Length() - i
		if n > 60 {
			n = 60
		}

		if len(highlights) > 0 {
			fmt.Fprintf(fp, "%16s", "")
			for _, h := range highlights[highlightIt:] {
				found := false
				if i+n >= h.Start && i < h.End {
					found = true
					for j := i; j < h.Start; j++ {
						fmt.Fprintf(fp, " ")
					}
					for j := h.Start; j < h.End && j < i+n; j++ {
						fmt.Fprintf(fp, "%c", h.Char)
					}
				}
				if i+n >= h.End {
					highlightIt++
				}
				if found { // Don't attempt to print overlapping highlights
					break
				}
			}
			fmt.Fprintf(fp, "\n")
		}

		for j, w := range which {
			nts := g.Nts[w][i : i+n]
			fmt.Fprintf(fp, "%-16s%s\t%d\n", names[j], string(nts), i+n)
		}

		if len(which) > 1 {
			fmt.Fprintf(fp, "%-16s%s\n", "", g.compare(i, i+n, which...))
		}

		if withTranslation {
			for j, _ := range which {
				var aas []byte
				nextIndex, aas = trans[j].TranslateLong(index, i, i+n)
				fmt.Fprintf(fp, "%-16s%s\n", names[j], string(aas))
			}
			fmt.Fprintf(fp, "\n")
		}
		index = nextIndex
	}

	fp.Flush()
	return nil
}

// Save in a clu-style format with the translation. Assume a and b are aligned.
// -1 for which means save everything.
func (g *Genomes) SaveWithTranslation(fname string,
	highlights []Highlight, which ...int) error {
	return g.saveCluStyle(fname, highlights, true, which...)
}

func (g *Genomes) SaveClu(fname string,
	highlights []Highlight, which ...int) error {
	return g.saveCluStyle(fname, highlights, false, which...)
}
