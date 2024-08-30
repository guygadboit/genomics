package genomes

import (
	"bufio"
	"errors"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"slices"
	"regexp"
	"strings"
)

/*
Represents a collection of aligned genomes (often two) with one genome per row.
The orfs "belong" to the first one in the set. We also use this for a single
genome, and sometimes for a protein
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

func (g *Genomes) CheckOrfs() error {
	var bad bool
	for _, orf := range g.Orfs {
		if orf.Start < 0 || orf.Start >= g.Length() {
			bad = true
			break
		}
		if orf.End < 0 || orf.End >= g.Length() {
			bad = true
			break
		}
	}
	if bad {
		return errors.New("Genome is too short for ORFs")
	}
	return nil
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

// Make a deep copy (of the nts-- the names are just shallow copied)
func (g *Genomes) Clone() *Genomes {
	ret := NewGenomes(g.Orfs, g.NumGenomes())
	for i := 0; i < len(g.Nts); i++ {
		ret.Nts[i] = make([]byte, len(g.Nts[i]))
		copy(ret.Nts[i], g.Nts[i])
	}
	for i := 0; i < g.NumGenomes(); i++ {
		ret.Names[i] = g.Names[i]
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

// Make a shallow copy of g but with a and b swapped. You could do this with
// filter but it would be annoying.
func (g *Genomes) Swap(a, b int) *Genomes {
	ret := NewGenomes(g.Orfs, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret.Nts[i] = g.Nts[i]
		ret.Names[i] = g.Names[i]
	}
	ret.Nts[a], ret.Nts[b] = ret.Nts[b], ret.Nts[a]
	ret.Names[a], ret.Names[b] = ret.Names[b], ret.Names[a]
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

func Min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

/*
How many differences are there looking at length bytes from pos in a and
comparing them to length bytes from pos+offset in b?
*/
func countDifferences(a []byte, b []byte, pos, offset, length int) int {
	var ret int

	aSpace := len(a) - (pos + length)
	bSpace := len(b) - (pos + offset + length)
	space := Min(aSpace, bSpace)
	if space < 0 {
		length += space
	}

	for i := 0; i < length; i++ {
		if a[pos+i] != b[pos+offset+i] {
			ret++
		}
	}

	return ret
}

// Insert '-' into other, or cut pieces out of other, to make it match ref
func align(ref []byte, other []byte) []byte {
	// Assuming things don't match up at pos, would it help to put a gap into a
	// here?
	needGap := func(pos int) bool {
		// We need a gap if things aren't already back on track a little down
		// the road.
		return countDifferences(other, ref, pos+20, 0, 20) > 2
	}

	// Would putting a gap this long in to a at pos be a good idea?
	tryGap := func(pos, length int, a, b []byte) bool {
		return countDifferences(a, b, pos+20, length, 20) <= 2
	}

	// Actually put a gap in
	insertGap := func(pos, length int) {
		other = append(other[:pos],
			append(make([]byte, length), other[pos:]...)...)
		for i := pos; i < pos+length; i++ {
			other[i] = '-'
		}
	}

	removeSection := func(pos, length int) {
		other = append(other[:pos], other[pos+length:]...)
	}

	for i := 0; ; i++ {
		if i >= len(ref) || i >= len(other) {
			break
		}
		if other[i] == ref[i] {
			continue
		}
		if !needGap(i) {
			continue
		}

		for j := 1; j < 1000; j++ {
			if tryGap(i, j, other, ref) {
				insertGap(i, j)
				i += j
				break
			}
			if tryGap(i, j, ref, other) {
				removeSection(i, j)
				i += j
				break
			}
		}
	}
	return other
}

func (g *Genomes) AlignCombine(other *Genomes) error {
	for i := 0; i < other.NumGenomes(); i++ {

		other.Nts[i] = align(g.Nts[0], other.Nts[i])
		n := len(other.Nts[i])

		// If other is too long at this point, just skip it. We'll see how
		// often this happens.
		if n > len(g.Nts[0]) {
			log.Printf("Alignment of %d failed.\n", i)
			continue
		} else {
			for j := 0; j < len(g.Nts[0])-n; j++ {
				other.Nts[i] = append(other.Nts[i], '-')
			}
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

func (g *Genomes) ToVector(which int) []float64 {
	ret := make([]float64, g.Length())
	for i := 0; i < g.Length(); i++ {
		switch g.Nts[which][i] {
		case 'A':
			ret[i] = 1.0
		case 'G':
			ret[i] = 2.0
		case 'T':
			ret[i] = 3.0
		case 'C':
			ret[i] = 4.0
		}
	}
	return ret
}

/*
If each genome was a vector, where -, A, G, T, C == 0, 1, 2, 4, 5 return the
position of the centroid
*/
func (g *Genomes) Centroid() []float64 {
	ret := make([]float64, g.Length())
	for i := 0; i < g.NumGenomes(); i++ {
		utils.VecAdd(ret, g.ToVector(i))
	}
	n := float64(g.NumGenomes())
	for i := 0; i < g.Length(); i++ {
		ret[i] /= n
	}
	return ret
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

func ParseHighlights(s string, sep string,
	oneBased bool, hChar byte) []Highlight {
	if s == "" {
		return nil
	}
	positions := utils.ParseInts(s, sep)
	ret := make([]Highlight, len(positions))
	for i, pos := range positions {
		if oneBased {
			pos -= 1
		}
		ret[i] = Highlight{pos, pos + 1, hChar}
	}
	return ret
}

func ParseHighlightFile(fname string,
	oneBased bool, hChar byte) ([]Highlight, error) {
	ret := make([]Highlight, 0)
	var err error
	re := regexp.MustCompile(`[GACT](\d+)[GACT]`)
	utils.Lines(fname, func(line string, lineErr error) bool {
		if lineErr != nil {
			err = lineErr
			return false
		}

		// Make this just work if the file contains mutations rather than
		// positions.
		groups := re.FindStringSubmatch(line)
		if groups != nil {
			line = groups[1]
		}

		pos := utils.Atoi(line)
		if oneBased {
			pos -= 1
		}
		ret = append(ret, Highlight{pos, pos+1, hChar})
		return true
	})
	return ret, err
}

// Returns the string to use for the highlights. If you wanted to optimize this
// you could track where you got to in highlights and resume searching from
// there the next time (as it is sorted).
func makeHighlightString(highlights []Highlight, start, length int) string {
	ret := make([]byte, length)
	// Start off with it being all spaces
	for i := 0; i < len(ret); i++ {
		ret[i] = ' '
	}

	// Now find any highlights that might be in there and put them in
positions:
	for i := 0; i < len(ret); i++ {
		for _, h := range highlights {
			if start+i >= h.Start && start+i < h.End {
				ret[i] = h.Char
				// If more than one highlight matches this position, just stick
				// with the first one we found.
				continue positions
			}
		}
	}

	return string(ret)
}

func (g *Genomes) saveCluStyle(fname string,
	highlights []Highlight, withTranslation bool, which ...int) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()
	fp := bufio.NewWriter(fd)

	if len(which) == 0 {
		which = make([]int, g.NumGenomes())
		for i := 0; i < g.NumGenomes(); i++ {
			which[i] = i
		}
	}

	names := make([]string, len(which))
	trans := make([]Translation, len(which))
	var index, nextIndex int
	for i, w := range which {
		names[i] = shorten(g.Names[w], 16)
		if withTranslation {
			trans[i] = Translate(g, w)
		}
	}

	slices.SortFunc(highlights, func(a, b Highlight) int {
		return a.Start - b.Start
	})

	for i := 0; i < g.Length(); i += 60 {
		n := g.Length() - i
		if n > 60 {
			n = 60
		}

		if len(highlights) > 0 {
			fmt.Fprintf(fp, "%16s%s\n",
				"", makeHighlightString(highlights, i, n))
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
		}
		fmt.Fprintf(fp, "\n")
		index = nextIndex
	}

	fp.Flush()
	return nil
}

// Save in a clu-style format with the translation. Nothing for which means
// save everything
func (g *Genomes) SaveWithTranslation(fname string,
	highlights []Highlight, which ...int) error {
	return g.saveCluStyle(fname, highlights, true, which...)
}

func (g *Genomes) SaveClu(fname string,
	highlights []Highlight, which ...int) error {
	return g.saveCluStyle(fname, highlights, false, which...)
}

// Given a position (which skips gaps) into the "which"th genome, convert that
// into a position in the 0th. So basically put the gaps back.
func (g *Genomes) ConvertPosition(which int, pos int) int {
	var gaps int
	for i := 0; i < pos; i++ {
		if g.Nts[which][i] == '-' {
			fmt.Printf("gap at %d\n", i)
			gaps++
		}
	}
	return pos + gaps
}
