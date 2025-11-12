package genomes

import (
	"bufio"
	"errors"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"reflect"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

var CodonTable = map[string]byte{
	"TTT": 'F', // Phenylalanine
	"TTC": 'F',

	"TTA": 'L', // Leucine
	"TTG": 'L',
	"CTT": 'L',
	"CTC": 'L',
	"CTA": 'L',
	"CTG": 'L',

	"ATT": 'I', // Isoleucine
	"ATC": 'I',
	"ATA": 'I',

	"ATG": 'M', // Methionine

	"GTT": 'V', // Valine
	"GTC": 'V',
	"GTA": 'V',
	"GTG": 'V',

	"TCT": 'S', // Serine
	"TCC": 'S',
	"TCA": 'S',
	"TCG": 'S',

	"CCT": 'P', // Proline
	"CCC": 'P',
	"CCA": 'P',
	"CCG": 'P',

	"ACT": 'T', // Threonine
	"ACC": 'T',
	"ACA": 'T',
	"ACG": 'T',

	"GCT": 'A', // Alanine
	"GCC": 'A',
	"GCA": 'A',
	"GCG": 'A',

	"TAT": 'Y', // Tyrosine
	"TAC": 'Y',

	"TAA": '*', // Stop
	"TAG": '*',

	"CAT": 'H', // Histidine
	"CAC": 'H',

	"CAA": 'Q', // Glutadine
	"CAG": 'Q',

	"AAT": 'N', // Asparagine
	"AAC": 'N',

	"AAA": 'K', // Lysine
	"AAG": 'K',

	"GAT": 'D', // Aspartic acid
	"GAC": 'D',

	"GAA": 'E', // Glutamic acid
	"GAG": 'E',

	"TGT": 'C', // Cysteine
	"TGC": 'C',

	"TGA": '*', // Stop
	"TGG": 'W', // Tryptophan

	"CGT": 'R', // Arginine
	"CGC": 'R',
	"CGA": 'R',
	"CGG": 'R',

	"AGT": 'S', // Serine
	"AGC": 'S',

	"AGA": 'R', // Arginine (again)
	"AGG": 'R',

	"GGT": 'G', // Glycine
	"GGC": 'G',
	"GGA": 'G',
	"GGG": 'G',
}

var ReverseCodonTable map[byte][]string

type Orf struct {
	Start, End int
	Name       string
	Reverse    bool // Reverse Complement
}

type Orfs []Orf

func LoadOrfs(fname string) Orfs {
	ret := make(Orfs, 0)

	fd, err := os.Open(fname)
	if err != nil {
		log.Fatalf("Can't open file %s", fname)
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)

loop:
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

		var reverse bool
		pat := regexp.MustCompile(`complement\((.*)\)`)

		line = strings.TrimSpace(line)
		m := pat.FindSubmatch([]byte(line))
		if len(m) != 0 {
			reverse = true
			line = string(m[1])
		}

		fields := strings.Fields(line)

		if strings.Contains(fields[0], ".") {
			subFields := strings.FieldsFunc(fields[0], func(r rune) bool {
				return r == '.'
			})
			fields = append(subFields, fields[1:]...)
		}

		start, err := strconv.Atoi(fields[0])
		if err != nil {
			log.Fatalf("Parse error in ORFs: <%s>\n", line)
		}

		// ORFs seem to be conventionally 1-based
		start -= 1

		end, err := strconv.Atoi(fields[1])
		if err != nil {
			log.Fatalf("Parse error in ORFs: <%s>", line)
		}

		var name string
		if len(fields) == 3 {
			name = fields[2]
		}

		ret = append(ret, Orf{start, end, name, reverse})
	}

	return ret
}

/*
Return the start of the codon and where pos is in it. We just do a linear
search since there aren't usually that many ORFs and this is probably as
fast as anything else
*/
func (orfs Orfs) GetCodonOffset(pos int) (int, int, error) {
	for i := 0; i < len(orfs); i++ {
		orf := &orfs[i]
		if pos >= orf.Start && pos < orf.End {
			orfPos := pos - orf.Start // pos relative to start of ORF
			return orf.Start + (orfPos/3)*3, orfPos % 3, nil
		}
	}
	return 0, 0, errors.New("Not in ORF")
}

/*
Convert a position into an ORF relative one. Return the index of the ORF and
the position in it.
*/
func (orfs Orfs) GetOrfRelative(pos int) (int, int, error) {
	for i, orf := range orfs {
		if pos >= orf.Start && pos < orf.End {
			return i, pos - orf.Start, nil
		}
	}
	return 0, 0, nil
}

func (orfs Orfs) Find(name string) (Orf, error) {
	for _, o := range orfs {
		if o.Name == name {
			return o, nil
		}
	}
	return Orf{}, errors.New("Not Found")
}

// The "Environment" of a subsequence is the codon-aligned section that
// completely contains it.
type Environment struct {
	start  int // Index into the original genome
	length int // How many nts in the subsequence this represents

	window  []byte // The whole aligned section
	offset  int    // The offset to the start of the subsequence
	protein []byte // Its translation
}

// rounded up to the nearest multiple of 3
func ceil3(n int) int {
	return n + 3 - n%3
}

/*
Assume nts are codon aligned and return a translation, with one amino-acid
letter per nt, so something like LLLRRRIII
*/
func TranslateAligned(nts []byte) []byte {
	ret := make([]byte, len(nts))

	for i := 0; i < len(nts)-3+1; i += 3 {
		aa, there := CodonTable[string(nts[i:i+3])]
		if !there {
			aa = '-'
		}
		for j := 0; j < 3; j++ {
			ret[i+j] = aa
		}
	}
	return ret
}

/*
Assume nts are codon aligned and return a translation, with one amino-acid
letter per aa, so something like LRI
*/
func TranslateAlignedShort(nts []byte) []byte {
	ret := make([]byte, 0)

	for i := 0; i < len(nts)-3+1; i += 3 {
		aa, there := CodonTable[string(nts[i:i+3])]
		if !there {
			aa = '-'
		}
		ret = append(ret, aa)
	}
	return ret
}

func (env *Environment) Init(genome *Genomes,
	pos int, n int, which int) error {
	env.start = pos
	env.length = n

	windowStart, codonOffset, err := genome.Orfs.GetCodonOffset(pos)
	if err != nil {
		return err
	}

	windowLen := ceil3(codonOffset + n)
	windowEnd := windowStart + windowLen

	env.offset = codonOffset
	env.window = genome.Nts[which][windowStart:windowEnd]

	if !utils.IsRegularPattern(env.window) {
		// Usually because of a gap ('-') in an alignment
		return errors.New("Non-nt in sequence")
	}

	env.protein = TranslateAligned(env.window)
	return nil
}

/*
Actually rewrite the environment with replacement nts. Not to be confused with
"Replace" which tells you whether it would be silent if you did rewrite it (but
doesn't actually do so).
*/
func (env *Environment) Rewrite(replacement []byte) error {
	if len(replacement) != env.length {
		return errors.New("Replacement is the wrong length")
	}
	if !utils.IsRegularPattern(replacement) {
		return errors.New("Non-nt in replacement sequence")
	}

	newWindow := make([]byte, len(env.window))
	copy(newWindow, env.window)
	copy(newWindow[env.offset:env.offset+env.length], replacement)

	env.window = newWindow
	env.protein = TranslateAligned(env.window)
	return nil
}

func (env *Environment) Subsequence() []byte {
	return env.window[env.offset : env.offset+env.length]
}

func (env *Environment) Protein() []byte {
	return env.protein[env.offset : env.offset+env.length]
}

// The protein in more conventional format, showing each AA one at a time
// instead of in triplets.
func (env *Environment) ProteinShort() []byte {
	ret := make([]byte, 0)
	for i := env.offset; i < env.offset+env.length; i += 3 {
		ret = append(ret, env.protein[i])
	}
	return ret
}

func (env *Environment) Print() {
	fmt.Println(string(env.Subsequence()))
	fmt.Println(string(env.Protein()))
}

/*
If we were to replace the subsequence this is the environment of, would
that be silent, and how many mutations would it contain?
*/
func (env *Environment) Replace(replacement []byte) (bool, int) {
	altWindow := make([]byte, len(env.window))
	copy(altWindow, env.window)
	copy(altWindow[env.offset:env.offset+env.length], replacement)

	protein := env.protein
	altProtein := TranslateAligned(altWindow)

	silent := true
	for i := 0; i < len(protein); i += 3 {
		if altProtein[i] != protein[i] {
			silent = false
			break
		}
	}

	subseq := env.Subsequence()
	differences := 0
	for i := 0; i < env.length; i++ {
		if replacement[i] != subseq[i] {
			differences++
		}
	}

	return silent, differences
}

// Silent alternative to a sequence of nts and how many muts that would require
type Alternative struct {
	NumMuts int
	Nts     []byte

	// In the case where NumMuts is 1, are we changing the 1st, 2nd or 3rd nt
	// of a codon? We use 1,2,3 so that 0 can mean "don't know"
	CodonPos int
}

type Alternatives []Alternative

// Iterator for finding the alternatives to a given subsequence
type altIter struct {
	protein  []byte // the protein we're finding nts for, as RL not RRRLLL
	odometer []int  // tracks the codon combinations as we iterate them
}

func (it *altIter) Init(protein []byte) {
	it.protein = protein
	it.odometer = make([]int, len(protein))
}

/*
Return the next alternative and whether there are any more to come after
it.
*/
func (it *altIter) Next() ([]byte, bool) {
	prot := it.protein
	ret := make([]byte, 0, len(prot)*3)

	for i := 0; i < len(prot); i++ {
		codons := ReverseCodonTable[prot[i]]
		ret = append(ret, []byte(codons[it.odometer[i]])...)
	}

	// Increment the odometer like a sort of odometer
	for j := 0; j < len(prot); j++ {
		codons := ReverseCodonTable[prot[j]]
		if it.odometer[j]+1 < len(codons) {
			it.odometer[j]++
			for k := 0; k < j; k++ {
				it.odometer[k] = 0
			}
			return ret, true
		}
	}

	return ret, false
}

func TestAlternatives() {
	protein := []byte("RL")
	var it altIter
	it.Init(protein)

	for {
		alt, more := it.Next()
		fmt.Println(string(alt))
		if !more {
			break
		}
	}
}

/*
Newer versions of Go have a more "ergonomic" slices.SortFunc which saves
you doing all this.
*/
func (a Alternatives) Len() int {
	return len(a)
}

func (a Alternatives) Less(i, j int) bool {
	return a[i].NumMuts < a[j].NumMuts
}

func (a Alternatives) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

// Return a copy of a which has been uniqued on Nts
func (a Alternatives) Unique() Alternatives {
	m := make(map[string]Alternative)
	for _, alt := range a {
		m[string(alt.Nts)] = alt
	}

	ret := make(Alternatives, 0)
	for _, v := range m {
		ret = append(ret, v)
	}
	return ret
}

/*
Find the alternative nt sequences that would not change the protein here,
ordered by fewest muts first. If constrained, an alternative doesn't count if
it needs changes outside env to be silent.
*/
func (env *Environment) FindAlternatives(maxMuts int,
	constrained bool) Alternatives {
	var it altIter
	ret := make(Alternatives, 0)

	// The protein stored in env is like LLLRRRIII. We want just LRI.
	windowLen := len(env.window)
	protein := make([]byte, windowLen/3)
	for i := 0; i < len(protein); i++ {
		protein[i] = env.protein[i*3]
	}

	it.Init(protein)
	existing := env.Subsequence()

	for more := true; more; {
		var alt []byte
		alt, more = it.Next()
		start, end := env.offset, env.offset+env.length

		if constrained {
			if !reflect.DeepEqual(alt[:start], env.window[:start]) {
				continue
			}
			if !reflect.DeepEqual(alt[end:], env.window[end:]) {
				continue
			}
		}

		numMuts, codonPos := 0, 0
		for i := 0; i < env.length; i++ {
			if alt[start+i] != existing[i] {
				numMuts++
				codonPos = (env.offset+i)%3 + 1
			}
		}

		if numMuts != 1 {
			codonPos = 0 // meaningless in this case
		}

		if numMuts > 0 && numMuts <= maxMuts {
			ret = append(ret, Alternative{numMuts,
				alt[start:end], codonPos})
		}

		if !more {
			break
		}
	}

	sort.Sort(ret)
	return ret
}

/*
Just translate a whole genome, iterating over all the codons in it from the
start
*/
type CodonIter struct {
	genome *Genomes // Alignment of genomes
	which  int      // Which one you want to translate
	orfI   int      // Which ORF we're in
	pos    int      // Where we are in it
}

func (it *CodonIter) Init(genome *Genomes, which int) {
	it.genome = genome
	it.which = which
	it.pos = genome.Orfs[0].Start
}

func (it *CodonIter) Next() (pos int, codon string, aa byte, err error) {
	genome := it.genome
	nts := genome.Nts
	for ; it.orfI < len(genome.Orfs); it.orfI++ {
		orf := genome.Orfs[it.orfI]
		if it.pos == -1 {
			it.pos = orf.Start
		}
		pos = it.pos
		if pos+3 <= orf.End {
			if orf.Reverse {
				rpos := orf.End - (pos - orf.Start) - 3
				codon = string(utils.ReverseComplement(
					nts[it.which][rpos : rpos+3]))
			} else {
				codon = string(nts[it.which][pos : pos+3])
			}
			var there bool
			aa, there = CodonTable[codon]
			if !there {
				aa = '-'
			}
			err = nil
			it.pos = pos + 3
			return
		}
		it.pos = -1 // Marks that we just moved to a new ORF.
	}
	return 0, "", 0, errors.New("No more ORFs")
}

type Codon struct {
	Pos int
	Nts string
	Aa  byte
}

func (c *Codon) Init(pos int, nts string) {
	c.Pos = pos
	c.Nts = nts

	aa, there := CodonTable[c.Nts]
	if !there {
		aa = '-'
	}
	c.Aa = aa
}

func (c Codon) ToString() string {
	return fmt.Sprintf("%c (%s)", c.Aa, c.Nts)
}

type Translation []Codon

func Translate(genome *Genomes, which int) Translation {
	ret := make(Translation, 0)
	var ci CodonIter
	ci.Init(genome, which)

	for {
		pos, nts, aa, err := ci.Next()
		if err != nil {
			break
		}
		ret = append(ret, Codon{pos, nts, aa})
	}
	return ret
}

type TranslationMap map[int]Codon

func NewTranslationMap(trans Translation) TranslationMap {
	ret := make(TranslationMap)

	for _, codon := range trans {
		for i := 0; i < 3; i++ {
			ret[codon.Pos+i] = codon
		}
	}

	return ret
}

// Returns -1 for not there
func (t Translation) Find(protein []byte, start int) int {
	for i := start; i < len(t); i++ {
		found := true
		for j := 0; j < len(protein); j++ {
			if t[i+j].Aa != protein[j] {
				found = false
				break
			}
		}
		if found {
			return t[i].Pos
		}
	}
	return -1
}

/*
Return the translation between start and end in "long" form so RRR for an R
codon. Start looking from index, and return the codon index (use this to
improve performance if you're calling this in a loop for the whole genome).
It should always work just to pass 0 for index, but be slower.
*/
func (t Translation) TranslateLong(index, start, end int) (int, []byte) {
	ret := make([]byte, 0)

	for pos := start; pos < end; pos++ {
		for ; index < len(t) && t[index].Pos+2 < pos; index++ {
		}
		if index == len(t) {
			break
		}
		codon := t[index]

		offset := pos - codon.Pos
		if offset < 0 {
			ret = append(ret, '|')
			continue
		}
		if codon.Aa == 0 {
			ret = append(ret, 'X')
		} else {
			ret = append(ret, codon.Aa)
		}
	}

	return index, ret
}

func (t Translation) ToString() string {
	var ret string

	for _, c := range t {
		ret += string(c.Aa)
	}
	return ret
}

/*
Do genomes a and b in an alignment code for the same thing between pos and
pos+length? Return silent, and the number of muts
*/
func IsSilent(g *Genomes, pos int, length int, a, b int) (bool, int, error) {
	var envA, envB Environment

	numMuts := utils.NumMuts(g.Nts[a][pos:pos+length],
		g.Nts[b][pos:pos+length])

	err := envA.Init(g, pos, length, a)
	if err != nil {
		return false, numMuts, err
	}

	err = envB.Init(g, pos, length, b)
	if err != nil {
		return false, numMuts, err
	}

	silent := reflect.DeepEqual(envA.Protein(), envB.Protein())
	return silent, numMuts, nil
}

/*
If you replaced a with replacement at pos, would it be silent relative to b?
Returns silent and the number of muts
*/
func IsSilentWithReplacement(g *Genomes,
	pos int, a, b int, replacement []byte) (bool, int, error) {
	var envA, envB Environment
	length := len(replacement)

	numMuts := utils.NumMuts(g.Nts[a][pos:pos+length], replacement)

	err := envA.Init(g, pos, length, a)
	if err != nil {
		return false, numMuts, err
	}

	err = envA.Rewrite(replacement)
	if err != nil {
		return false, numMuts, err
	}

	err = envB.Init(g, pos, length, b)
	if err != nil {
		return false, numMuts, err
	}

	silent := reflect.DeepEqual(envA.Protein(), envB.Protein())
	return silent, numMuts, nil
}

// Returns whether this was silent and the old and new proteins at that location
func ProteinChange(g *Genomes,
	pos int, a, b int, replacement []byte) (bool, []byte, []byte, error) {
	var envA, envB Environment
	length := len(replacement)

	err := envA.Init(g, pos, length, a)
	if err != nil {
		return false, nil, nil, err
	}

	err = envA.Rewrite(replacement)
	if err != nil {
		return false, nil, nil, err
	}

	err = envB.Init(g, pos, length, b)
	if err != nil {
		return false, nil, nil, err
	}

	silent := reflect.DeepEqual(envA.Protein(), envB.Protein())
	return silent, envB.Protein(), envA.Protein(), nil
}

func newReverseTable() map[byte][]string {
	ret := make(map[byte][]string)

	for k := range CodonTable {
		v, _ := CodonTable[k]
		codons, there := ret[v]
		if !there {
			codons = make([]string, 0)
		}
		ret[v] = append(codons, k)
	}
	return ret
}

func init() {
	ReverseCodonTable = newReverseTable()
}
