package main

import (
	"bufio"
	"encoding/gob"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
	"io"
	"log"
	"os"
	"reflect"
	"regexp"
	"sort"
	"strings"
)

type Direction int

const (
	UNKNOWN Direction = iota
	FORWARDS
	BACKWARDS
)

type Homology struct {
	Forwards  int     // number of nts that match ahead of the insert
	Backwards int     // number of nts that match behind it
	E         float64 // E value of the match with the most homology
	FullMatch []byte  // The match including the homology
}

type Match struct {
	stats.BlastResult
	Name      string
	SourcePos int
	Dir       Direction
	Homology  Homology
	// XHomology  Homology // weird "crossed-over" homology
}

type Insertion struct {
	Id            int               // The order we found them in. Should be a line number
	Pos           utils.OneBasedPos // Where
	Nts           []byte            // What
	NSeqs         int               // How many times
	InWH1         bool              // Is it in WH1?
	InHuman       bool              // Is it in human? Note that most short things will be
	PosInHuman    int
	DirInHuman    Direction
	InOrf         bool // Is it in an ORF?
	NumHere       int  // Number of insertions at this location
	StrictNumHere int  // Number here with NSeqs > 1
	Matches       []Match
}

func (i *Insertion) ToString() string {
	return fmt.Sprintf("%s at %d (%d seqs)", string(i.Nts),
		i.Pos, i.NSeqs)
}

type InsertionData struct {
	Locations       map[utils.OneBasedPos]int // maps position to count
	LocationsStrict map[utils.OneBasedPos]int // require nseqs > 1
	Insertions      []Insertion
	NucDistro       *mutations.NucDistro
}

type InsertionNtIterator struct {
	id      *InsertionData
	filters []filterFunc
	i, j    int
}

func (it *InsertionNtIterator) Start() {
	it.i = 0
}

func (it *InsertionNtIterator) Get() byte {
	return it.id.Insertions[it.i].Nts[it.j]
}

func filtersPass(ins *Insertion, filters []filterFunc) bool {
	for _, f := range filters {
		if !f(ins) {
			return false
		}
	}
	return true
}

func (it *InsertionNtIterator) Next() {
	it.j++
	if it.j == len(it.id.Insertions[it.i].Nts) {
		it.j = 0

		for {
			it.i++
			if it.i == len(it.id.Insertions) {
				break
			}
			if filtersPass(&it.id.Insertions[it.i], it.filters) {
				break
			}
		}
	}
}

func (it *InsertionNtIterator) End() bool {
	return it.i == len(it.id.Insertions)
}

func (id *InsertionData) GetNucDistro(
	filters []filterFunc) *mutations.NucDistro {
	if id.NucDistro != nil {
		return id.NucDistro
	}
	it := InsertionNtIterator{id, filters, 0, 0}
	id.NucDistro = mutations.NewNucDistro(&it, mutations.NT_ALPHABET)
	fmt.Println(id.NucDistro)
	return id.NucDistro
}

// Keep the lengths and locations the same, randomize the actual nucleotides.
// Use a nucleotide distribution based on the insertions that match the
// filters.
func (id *InsertionData) Randomize(filters []filterFunc) {
	nd := id.GetNucDistro(filters)
	for i, _ := range id.Insertions {
		ins := &id.Insertions[i]
		nd.RandomSequence(ins.Nts)
	}
}

func LoadInsertions(fname string, minLen int, minSeqs int) []Insertion {
	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal("Can't open %s", fname)
	}
	defer fd.Close()

	ret := make([]Insertion, 0)
	fp := bufio.NewReader(fd)

	// Match for [GATC]+, so ignore any with Ns or weird ambiguous nts like HDK
	// etc.
	pat := regexp.MustCompile(`(\d+) ins_(\d+):([GATC]+) \((\d+) seqs\)`)

reading:
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break reading
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)
		groups := pat.FindAllStringSubmatch(line, -1)
		if groups == nil {
			continue
		}

		ins := Insertion{utils.Atoi(groups[0][1]),
			utils.OneBasedPos(utils.Atoi(groups[0][2])),
			[]byte(groups[0][3]),
			utils.Atoi(groups[0][4]),
			false, false,
			0, UNKNOWN, true, 0, 0, make([]Match, 0)}

		if len(ins.Nts) < minLen {
			continue
		}

		if ins.NSeqs < minSeqs {
			continue
		}

		// Some of the insertions seem to be nearly all A. This looks bogus. So
		// filter them out.
		var aCount int
		for i := 0; i < len(ins.Nts); i++ {
			if ins.Nts[i] == 'A' {
				aCount++
			}
		}
		if (float64(aCount) / float64(len(ins.Nts))) > 0.99 {
			continue
		}

		ret = append(ret, ins)
	}

	return ret
}

func (id *InsertionData) Find(fname string, minLen int, minSeqs int) {
	id.Insertions = LoadInsertions(fname, minLen, minSeqs)
	id.Locations = byLocation(id.Insertions, 1)
	id.LocationsStrict = byLocation(id.Insertions, 2)

	for i, v := range id.Insertions {
		id.Insertions[i].NumHere = id.Locations[v.Pos]
		id.Insertions[i].StrictNumHere = id.LocationsStrict[v.Pos]
	}
}

func (id *InsertionData) Save(fname string) {
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	enc := gob.NewEncoder(fp)
	err = enc.Encode(id)
	if err != nil {
		log.Fatal(err)
	}
	fp.Flush()
	fmt.Printf("Saved to %s\n", fname)
}

func (id *InsertionData) Load(fname string) {
	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)
	dec := gob.NewDecoder(fp)
	err = dec.Decode(id)
	if err != nil {
		log.Fatal(err)
	}
}

func (id *InsertionData) OutputInsertions(fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer f.Close()
	w := bufio.NewWriter(f)

	fmt.Fprintf(w, "id pos pattern nseqs"+
		" in_wh1 in_human num_here strict_num_here num_matches\n")
	fmt.Fprintf(w, "int int str int bool bool int int int\n")

	for _, ins := range id.Insertions {
		fmt.Fprintf(w, "%d %d %s %d %t %t %d %d %d\n", ins.Id,
			ins.Pos, string(ins.Nts), ins.NSeqs, ins.InWH1,
			ins.InHuman, ins.NumHere, ins.StrictNumHere, len(ins.Matches))
	}

	w.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func OutputMatchHeader(w *bufio.Writer) {
	fmt.Fprintf(w, "id name pattern forwards full_match seqs "+
		"num_here strict_num_here pos src_pos in_human in_wh1 forwards_h "+
		"backwards_h score E hE\n")
	fmt.Fprintf(w, "int str str bool str int int int int int "+
		"bool bool int int float float float\n")
}

func (m *Match) Output(w *bufio.Writer, ins *Insertion) {
	fullMatch := "-"
	if m.Homology.FullMatch != nil {
		fullMatch = string(m.Homology.FullMatch)
	}

	forwards := m.Dir == FORWARDS
	fmt.Fprintf(w, "%d %s %s %t %s %d %d %d %d %d "+
		"%t %t %d %d %f %g %g\n",
		ins.Id, m.Name, string(ins.Nts), forwards, fullMatch,
		ins.NSeqs, ins.NumHere, ins.StrictNumHere,
		ins.Pos, m.SourcePos, ins.InHuman,
		ins.InWH1, m.Homology.Forwards, m.Homology.Backwards,
		m.Score, m.E, m.Homology.E)
}

func OutputMatches(insertions []Insertion, fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer f.Close()
	w := bufio.NewWriter(f)

	OutputMatchHeader(w)
	for _, ins := range insertions {
		for _, m := range ins.Matches {
			m.Output(w, &ins)
		}
	}
	w.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func (id *InsertionData) OutputMatches(fname string) {
	OutputMatches(id.Insertions, fname)
}

func Summary(insertions []Insertion) {
	fmt.Printf("%d insertions\n", len(insertions))

	var total int
	for i := 0; i < len(insertions); i++ {
		total += len(insertions[i].Nts)
	}
	fmt.Printf("%d nts altogether (average length %.2f)\n",
		total, float64(total)/float64(len(insertions)))
}

/*
Call cb for insertion if it's found in nts forwards or backwards
*/
func Search(ins *Insertion, index string, cb func(*Insertion, int, bool)) {
	for s := genomes.NewBidiIndexSearch(index, ins.Nts); !s.End(); s.Next() {
		pos, _ := s.Get()
		cb(ins, pos, s.IsForwards())
	}
}

func findInVirus(insertions []Insertion, minLength, maxLength int) {
	reportFound := func(ins *Insertion) {
		ins.InWH1 = true
	}

	var found, count int

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]

		nts := ins.Nts
		if len(nts) < minLength {
			continue
		}
		if len(nts) >= maxLength {
			continue
		}
		count++

		Search(ins, "/fs/f/genomes/viruses/SARS2/index", func(ins *Insertion,
			pos int, forwards bool) {
			reportFound(ins)
			found++
		})
	}

	fmt.Printf("Min length %d: %d (/%d) were found in SARS-CoV-2\n", minLength,
		found, count)
}

func loadHuman() *genomes.Genomes {
	fmt.Printf("Loading...\n")
	g := genomes.LoadGenomes("/fs/f/genomes/human/GRCh38_latest_genomic.fna.gz",
		"", true)
	fmt.Printf("Loaded human\n")
	return g
}

/*
Mark the ones you find in human, only if they aren't already in WH1
*/
func findInHuman(insertions []Insertion, minLength, maxLength int) {
	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		nts := ins.Nts

		if len(nts) < minLength {
			continue
		}

		if len(nts) >= maxLength {
			continue
		}

		if ins.InWH1 {
			continue
		}

		Search(ins, "/fs/f/genomes/human/index", func(ins *Insertion,
			pos int, forwards bool) {
			fmt.Printf("%d (length %d) is in human\n", ins.Id, len(ins.Nts))
			ins.InHuman = true
			if forwards {
				ins.DirInHuman = FORWARDS
			} else {
				ins.DirInHuman = BACKWARDS
			}
			// FIXME you could end the search as soon as you find one. But we
			// won't do this now as it's running and has already got quite
			// far... And it doesn't cost much.
		})

		if i%10 == 0 {
			fmt.Printf("Searched %d/%d\n", i, len(insertions))
		}
	}
}

func makeIndex(insertions []Insertion) map[int]*Insertion {
	insMap := make(map[int]*Insertion, len(insertions))

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		insMap[ins.Id] = ins
	}

	return insMap
}

/*
Marking those that are in human or not is very time-consuming, so reload
our saved results from a file. And also load which ones are in orfs while
we're at it.
*/
func loadInHuman(insertions []Insertion, insMap map[int]*Insertion) {
	fp := utils.NewFileReader("human.txt")
	defer fp.Close()

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

		line = strings.TrimSpace(line)
		fields := strings.Fields(line)

		id := utils.Atoi(fields[0])
		ins := insMap[id]

		if ins != nil && fields[2] == "human" {
			ins.InHuman = true
			ins.PosInHuman = utils.Atoi(fields[3])

			switch fields[4] {
			case "forwards":
				ins.DirInHuman = FORWARDS
			case "backwards":
				ins.DirInHuman = BACKWARDS
			}
		}
	}
}

/*
	Ponderous sort but less trouble than updating my whole system to the latest
	Go version that has slices.SortFunc
*/

type posCount struct {
	pos   utils.OneBasedPos
	count int
}

type posCounts []posCount

func (p posCounts) Len() int {
	return len(p)
}

func (p posCounts) Less(i, j int) bool {
	return p[i].pos < p[j].pos
}

func (p posCounts) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func byLocation(insertions []Insertion, minSeqs int) map[utils.OneBasedPos]int {
	positions := make(map[utils.OneBasedPos]int)

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		if ins.NSeqs >= minSeqs {
			count, _ := positions[ins.Pos]
			positions[ins.Pos] = count + 1
		}
	}
	return positions
}

func SaveLocations(insertions []Insertion, minLength int) map[utils.OneBasedPos]int {
	positions := byLocation(insertions, 0)

	posCounts := make(posCounts, 0, len(positions))
	for k, v := range positions {
		posCounts = append(posCounts, posCount{k, v})
	}

	sort.Sort(posCounts)

	fd, err := os.Create("locations.txt")
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	for i := 0; i < posCounts.Len(); i++ {
		pc := posCounts[i]
		fmt.Fprintf(w, "%d %d\n", pc.pos, pc.count)
	}

	w.Flush()
	fmt.Printf("Wrote locations.txt\n")
	return positions
}

type filterFlag int

const (
	ANYTHING      filterFlag = 0
	EXCLUDE_WH1              = 1
	EXCLUDE_HUMAN            = 1 << 1
)

// Return true if we are including this insertion, i.e. not filtering it out
type filterFunc func(*Insertion) bool

func makeMinLengthFilter(minLength int) filterFunc {
	return func(ins *Insertion) bool {
		return len(ins.Nts) >= minLength
	}
}

func makeMaxLengthFilter(maxLength int) filterFunc {
	return func(ins *Insertion) bool {
		return len(ins.Nts) <= maxLength
	}
}

func makeMinNumHereFilter(minimum int) filterFunc {
	return func(ins *Insertion) bool {
		return ins.NumHere >= minimum
	}
}

func makeMinStrictNumHereFilter(minimum int) filterFunc {
	return func(ins *Insertion) bool {
		return ins.StrictNumHere >= minimum
	}
}

func makePositionFilter(minPos, maxPos utils.OneBasedPos) filterFunc {
	return func(ins *Insertion) bool {
		if ins.Pos <= minPos {
			return false
		}
		if ins.Pos >= maxPos {
			return false
		}
		return true
	}
}

func makeFlagFilter(flags filterFlag) filterFunc {
	return func(ins *Insertion) bool {
		if flags&EXCLUDE_WH1 != 0 && ins.InWH1 {
			return false
		}
		if flags&EXCLUDE_HUMAN != 0 && ins.InHuman {
			return false
		}
		return true
	}
}

func makeMinSeqsFilter(minSeqs int) filterFunc {
	return func(ins *Insertion) bool {
		return ins.NSeqs >= minSeqs
	}
}

func makeCodonAlignFilter() filterFunc {
	return func(ins *Insertion) bool {
		return len(ins.Nts)%3 == 0
	}
}

func makeNotInWH1Filter() filterFunc {
	// This is assuming your insertions aren't already marked with inWH1. If
	// they are you can use the flag filter for more speed.
	return func(ins *Insertion) bool {
		for s := genomes.NewBidiIndexSearch("/fs/f/genomes/viruses/SARS2/index",
			ins.Nts); !s.End(); s.Next() {
			return false
		}
		return true
	}
}

// Does s contain count or more repeats of length length?
func HasRepeats(s []byte, count, length int) bool {
	for offset := 0; offset < length; offset++ {
		var (
			prev []byte
			n    int
		)
		for i := offset; i <= len(s)-length-offset; i += length {
			sl := s[i : i+length]
			if prev != nil {
				if reflect.DeepEqual(prev, sl) {
					n++
				} else {
					n = 0
				}
			}
			if n == count-1 {
				return true
			}
			prev = sl
		}
	}
	return false
}

func makeSillyFilter() filterFunc {
	pat := regexp.MustCompile(`G{6}|A{6}|T{6}|C{6}`)

	return func(ins *Insertion) bool {
		if pat.Find(ins.Nts) != nil {
			return false
		}

		if HasRepeats(ins.Nts, 5, 2) {
			return false
		}

		return true
	}
}

/*
Call cb on all the insertions that match the filter. This requires that
you've already marked which ones are in WH1 with findInVirus.
*/
func filterInsertions(insertions []Insertion,
	filters []filterFunc, cb func(*Insertion), verbose bool) {
	var name string

insertions:
	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]

		if filters != nil {
			for j := 0; j < len(filters); j++ {
				if !filters[j](ins) {
					continue insertions
				}
			}
		}

		if verbose {
			fmt.Printf("%s: %d ins_%d %s (%d seqs)\n", name,
				len(ins.Nts), ins.Pos, string(ins.Nts), ins.NSeqs)
		}

		cb(ins)
	}
}

/*
Output them all in one fasta file with NNN between each of them. Useful if
you want to look at dinucleotide composition etc.
*/
func OutputCombinedFasta(fname string, name string,
	insertions []Insertion,
	filters []filterFunc,
	verbose bool) int {
	sep := []byte("NNN")
	nts := make([]byte, 0)
	var count int

	cb := func(ins *Insertion) {

		nts = append(nts, ins.Nts...)
		nts = append(nts, sep...)
		count++
	}

	filterInsertions(insertions, filters, cb, verbose)

	var orfs genomes.Orfs
	genomes := genomes.NewGenomes(orfs, 1)
	genomes.Nts[0] = nts
	genomes.Names[0] = name
	genomes.Save("SC2 Insertions", fname, 0)

	fmt.Printf("%s: %d insertions\n", name, count)
	return count
}

/*
Put all the insertions matching the filters into a single fasta file with
headers between each one. I think you should then be able to BLAST the
whole lot.
*/
func OutputFasta(fname string, insertions []Insertion,
	filters []filterFunc, verbose bool) {

	var orfs genomes.Orfs
	genomes := genomes.NewGenomes(orfs, 0)

	cb := func(ins *Insertion) {
		genomes.Nts = append(genomes.Nts, ins.Nts)
		name := fmt.Sprintf("ins_%d_%d", ins.Pos, ins.Id)
		genomes.Names = append(genomes.Names, name)
	}

	filterInsertions(insertions, filters, cb, verbose)
	genomes.SaveMulti(fname)
}

func outputDinucs(fname string,
	insertions []Insertion,
	filters []filterFunc,
	verbose bool) {

	fd, err := os.Create(fname)
	if err != nil {
		log.Fatalf("Can't create %s", fname)
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)
	fmt.Fprintf(w, "id len G+C, CpG TpA CpGF insert human\n")

	cb := func(ins *Insertion) {
		dp := genomes.CalcProfile(ins.Nts)

		human := "unknown-if-human"
		if len(ins.Nts) >= 20 {
			if ins.InHuman {
				human = "human"
			} else {
				human = "not-human"
			}
		}

		fmt.Fprintf(w, "%d %d %.3f %.3f %.3f %.3f %s %s\n",
			ins.Id, len(ins.Nts), dp.GC, dp.CpG, dp.TpA, dp.CpGF,
			string(ins.Nts), human)
	}
	filterInsertions(insertions, filters, cb, verbose)
	w.Flush()
}

func showLength(insertions []Insertion) {
	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		if ins.InWH1 {
			continue
		}
		fmt.Printf("%d ins_%d %s\n", len(ins.Nts), ins.Pos, ins.Nts)
	}
}

/*
func otherHCoVs(insertions []Insertion, tol float64) {
	covs := [...]string{"229E", "NL63", "OC43", "HKU1", "WH1"}

	var wg sync.WaitGroup
	for i := 0; i < len(covs); i++ {
		wg.Add(1)
		go func(i int) {
			findInVirus(covs[i], insertions, alse, tol)
			wg.Done()
		}(i)
	}
	wg.Wait()
}
*/

func showHuman(insertions []Insertion) {
	g := loadHuman()
	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]

		if !ins.InHuman {
			continue
		}

		pos := ins.PosInHuman

		var dir string
		switch ins.DirInHuman {
		case FORWARDS:
			dir = "forwards"
		case BACKWARDS:
			dir = "backwards"
		}

		fmt.Printf("%d ins_%d (%d nts %d sequences) found at %d %s\n",
			ins.Id, ins.Pos, len(ins.Nts), ins.NSeqs, pos, dir)

		fmt.Println(string(ins.Nts))

		var comparison []byte
		if ins.DirInHuman == BACKWARDS {
			rc := utils.ReverseComplement(ins.Nts)
			fmt.Println(string(rc))
			comparison = rc
		} else {
			comparison = ins.Nts
		}

		humanBit := g.Nts[0][pos : pos+len(ins.Nts)]

		for i := 0; i < len(comparison); i++ {
			if comparison[i] == humanBit[i] {
				fmt.Printf("%c", '|')
			} else {
				fmt.Printf("%c", ' ')
			}
		}
		fmt.Printf("\n")
		fmt.Println(string(humanBit))
	}
}

func getCpG(ins *Insertion) float64 {
	dp := genomes.CalcProfile(ins.Nts)
	return dp.CpG
}

func appendFCS(insertions []Insertion) []Insertion {
	lastId := 0

	if len(insertions) > 0 {
		lastIns := &insertions[len(insertions)-1]
		lastId = lastIns.Id
	}

	insertions = append(insertions,
		Insertion{lastId + 1, 23601, []byte("CTCCTCGGCGGG"),
			2, false, false, 0, UNKNOWN, false, 50, 50, nil})

	return insertions
}

/*
Write in-orf.txt which marks which insertions are in ORFs. We will load
that back in later. It's a bit slow to keep doing it every time.
*/
func findOrfs(insertions []Insertion) {
	orfs := genomes.LoadOrfs("../fasta/WH1.orfs")
	fd, _ := os.Create("in-orf.txt")
	defer fd.Close()

	w := bufio.NewWriter(fd)

	for i := 0; i < len(insertions); i++ {
		ins := &insertions[i]
		_, _, err := orfs.GetCodonOffset(int(ins.Pos) - 1)
		in_orf := err == nil
		fmt.Fprintf(w, "%d %t\n", ins.Id, in_orf)
	}

	w.Flush()
}

func loadInOrfs(insertions []Insertion, insMap map[int]*Insertion) {
	fp := utils.NewFileReader("in-orf.txt")
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

		line = strings.TrimSpace(line)
		fields := strings.Fields(line)

		id := utils.Atoi(fields[0])
		ins := insMap[id]

		if ins != nil {
			insMap[id].InOrf = fields[1] == "true"
		}
	}
}

func testBlast() {
	bc := stats.BlastDefaultConfig()
	results := stats.Blast(bc, "bacteria/Treponema",
		[]byte("GCGGTGGAGCATGGGGTTTAATTCG"), 1e-4, 1, stats.VERBOSE)
	fmt.Println(results)
}

// Return the number of occurrences of alternative encodings of the FCS per
// million nucleotides.
func CountFCSAlternatives(wh1 *genomes.Genomes, source *Source) float64 {
	if wh1 == nil {
		wh1 = genomes.LoadGenomes("../fasta/WH1.fasta",
			"../fasta/WH1.orfs", false)
	}

	// Count the actual FCS first, just to be sure
	count, n := 0, 0
	search := genomes.NewBidiIndexSearch(
		source.Index, []byte("CTCCTCGGCGGG"))
	n += search.GenomeLength()
	for ; !search.End(); search.Next() {
		count++
	}
	fmt.Printf("%s actual: %d\n", source.Name, count)
	actual := float64(count*1e6) / float64(n)

	// Now count the alternative encodings
	var env genomes.Environment
	env.Init(wh1, 23600, 12, 0)
	alternatives := env.FindAlternatives(12, true)

	total, n := 0, 0
	for _, alt := range alternatives {
		count := 0
		search := genomes.NewBidiIndexSearch(
			source.Index, alt.Nts)
		n += search.GenomeLength()
		for count = 0; !search.End(); search.Next() {
			count++
		}
		fmt.Printf("%s %s: %d\n", source.Name, string(alt.Nts), count)
		total += count
	}

	ret := float64(total*1e6) / float64(n)
	fmt.Printf("%s: actual=%f alts=%f\n", source.Name, actual, ret)
	return ret
}

func CountPattern(sources []Source, pattern []byte) {
	fname := fmt.Sprintf("%s.txt", string(pattern))
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)

	for _, source := range sources {
		search := genomes.NewBidiIndexSearch(source.Index, pattern)
		var count int
		for ; !search.End(); search.Next() {
			count++
		}

		n := search.GenomeLength()
		freq := float64(count*1e6) / float64(n)
		fmt.Printf("%s: %d/%d %g per million\n",
			source.Name, count, n, freq)
		fmt.Fprintf(w, "%s %g\n", source.Name, freq)
	}

	w.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func CountAndSave(id *InsertionData) {
	id.Find("insertions2.txt", 6, 1)

	insertions := id.Insertions
	utils.Sort(len(insertions), true,
		func(i, j int) {
			insertions[i], insertions[j] = insertions[j], insertions[i]
		},
		nil,
		func(i int) float64 {
			return float64(len(insertions[i].Nts))
		})

	findInVirus(id.Insertions, 12, 200)

	/*
		filters := []filterFunc{
			makeMinLengthFilter(12),
			makeMaxLengthFilter(24),
		}

		fmt.Printf("Matching insertions\n")
		sources := GetSources(BACTERIA)
		CountInGenomes(nil, insertions, sources, filters, 0.0, BLAST|APPEND)
	*/
	// id.Save("insertions2.gob")
}

// Note that if randomize it randomizes the actual insertion data
func AddSource(data *InsertionData,
	name string, save bool, iterations int,
	tol float64, actions MatchAction) {
	sources := make([]Source, 0)

	all := GetSources(ANIMAL | BACTERIA)
	for _, source := range all {
		if source.Name == name {
			sources = append(sources, source)
		}
	}

	if len(sources) == 0 {
		log.Fatal("Can't find %s\n", name)
	}

	filters := []filterFunc{
		makeMinLengthFilter(12),
		makeMaxLengthFilter(24),
		makePositionFilter(0, 29870),
		makeSillyFilter(),
		makeCodonAlignFilter(),
		makeFlagFilter(EXCLUDE_WH1),
	}

	CountInGenomes(nil, data, sources, filters, tol, iterations, actions)

	if save {
		data.Save("insertions-new.gob")
		fmt.Printf("Wrote insertions-new.gob\n")
	}
}

func BlastFCS(sources []Source) {
	pattern := []byte("TTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAA")
	fcs := 20
	fcsLen := 12
	bc := stats.BlastDefaultConfig()
	for leftMargin := 0; leftMargin < 20; leftMargin++ {
		for rightMargin := 0; rightMargin < 20; rightMargin++ {
			blastPat := pattern[fcs-leftMargin : fcs+fcsLen+rightMargin]
			for _, source := range sources {
				results := stats.Blast(bc, source.Path,
					blastPat, 1, 1, stats.NOT_VERBOSE)
				if len(results) == 1 {
					fmt.Printf("%s %d %s %d len=%d %g\n", source.Name,
						leftMargin,
						string(blastPat),
						rightMargin,
						results[0].Length, results[0].E)
				}
			}
		}
	}
}

func UpdateHomology(wh1 *genomes.Genomes,
	insertions []Insertion, sources []Source) {
	fmt.Println("Updating homology...")
	if wh1 == nil {
		wh1 = genomes.LoadGenomes("../fasta/WH1.fasta",
			"../fasta/WH1.orfs", false)
	}

	sourceMap := make(map[string]*Source)
	for i, _ := range sources {
		source := &sources[i]
		sourceMap[source.Name] = source
	}

	for i, _ := range insertions {
		ins := &insertions[i]
		/*
			if ins.Id != 12244 {
				continue
			}
		*/
		for _, match := range ins.Matches {
			/*
				if match.Name != "ANaesl" || match.SourcePos != 3055679 {
					continue
				}
			*/
			source := sourceMap[match.Name]

			match.Homology = FindHomology(wh1, ins,
				source, match.SourcePos, match.Dir == FORWARDS, true)
		}
		if i%100 == 0 {
			fmt.Printf("Done %d/%d\n", i, len(insertions))
		}
	}
}

func findExpectedHomology() {
	root := "/fs/f/genomes/bacteria/"
	/*
		bacteria := make([]string, 3)
		for i, s := range []string{"ANaesl", "AVisc", "AIsrael"} {
			bacteria[i] = root + s + fmt.Sprintf("/%s.fasta.gz", s)
		}

		fmt.Println("All Actinomyces")
		FindExpectedHomologyFreq(int(1e6), bacteria...)

		for _, b := range bacteria {
			fmt.Println(b)
			FindExpectedHomologyFreq(int(3e6), b)
		}
	*/

	root = "/fs/f/genomes/"
	animals := []string{"cod", "human", "pangolin", "rabbit", "bat", "lizard"}
	for _, animal := range animals {
		fmt.Println(animal)
		FindExpectedHomologyFreq(int(3e6),
			root+animal+fmt.Sprintf("/%s.fasta.gz", animal))
	}
}

func ShowCGGInsertions(id *InsertionData) {
	for _, ins := range id.Insertions {
		if strings.Contains(string(ins.Nts), "CGGCGG") {
			fmt.Printf("%d %s\n", len(ins.Nts), ins.ToString())
		}
	}
}

func main() {
	var (
		data                 InsertionData
		countCGG, blastFCS   bool
		outputName           string
		expectedHomology     bool
		findHomology         bool
		fcsAlternatives      bool
		countFCSAlternatives bool
		outputGob            bool
		blast                bool
		blastHomology        bool
		randomize            bool
		iterations           int
		outputFasta          bool
		tol					float64
	)

	flag.BoolVar(&countCGG, "cgg", false, "Count CGGCGG")
	flag.BoolVar(&blastFCS, "fcs", false, "Blast FCS with context")
	flag.BoolVar(&expectedHomology, "calc-ef", false, "Calculate expected"+
		"homology frequency")
	flag.StringVar(&outputName, "output", "", "Output")
	flag.BoolVar(&findHomology, "homol", false, "Update homology")
	flag.BoolVar(&fcsAlternatives,
		"fcs-alt", false, "Look at FCS alternatives")
	flag.BoolVar(&countFCSAlternatives,
		"count-fcs-alt", false, "Count FCS alternatives")
	flag.BoolVar(&outputGob, "output-gob", false, "Output the gob")
	flag.BoolVar(&blast, "blast", false, "Blast the matches")
	flag.BoolVar(&blastHomology,
		"blast-h", false, "Blast the matches with homology")
	flag.BoolVar(&randomize,
		"random", false, "Randomize the insertions (!)")
	flag.IntVar(&iterations,
		"iterations", 1, "If randomizing how many times to do it")
	flag.BoolVar(&outputFasta,
		"output-fasta", false, "Output a fasta file")
	flag.Float64Var(&tol, "tol", 0.0, "Tolerance")
	flag.Parse()

	if _, err := os.Stat("insertions2.gob"); err == nil {
		data.Load("insertions2.gob")
	} else {
		CountAndSave(&data)
	}

	if outputName != "" {
		var actions MatchAction = OUTPUT

		if blast {
			actions |= BLAST
		}
		if blastHomology {
			actions |= BLAST_HOMOLOGY
		}
		if randomize {
			actions |= RANDOMIZE
		} else if iterations > 1 {
			log.Fatal("More than one iteration with randomize is pointless")
		}
		AddSource(&data, outputName, false, iterations, tol, actions)
	}

	if blastFCS {
		sources := GetSources(BACTERIA)
		BlastFCS(sources)
	}

	if expectedHomology {
		findExpectedHomology()
		return
	}

    if outputFasta {
        filters := []filterFunc{
            makeMinLengthFilter(12),
            makeMaxLengthFilter(24),
            makeMinSeqsFilter(1),
            makeMinStrictNumHereFilter(1),
            makeSillyFilter(),
            makeCodonAlignFilter(),
            makePositionFilter(0, 29870),
        }
        OutputFasta("../fasta/SplitInsertions.fasta",
			data.Insertions, filters, false)
        OutputCombinedFasta("../fasta/CombinedInsertions.fasta",
			"Insertions", data.Insertions, filters, false)
    }

	/*
		findInVirus(data.Insertions, 12, 200)
		findInHuman(data.Insertions, 12, 200)
		data.Save("insertions2b.gob")
	*/

	if findHomology {
		sources := GetSources(BACTERIA)
		UpdateHomology(nil, data.Insertions, sources)
		data.Save("insertions-new.gob")
		fmt.Println("Wrote insertions-new.gob with updated homologies")
	}

	if countCGG {
		ShowCGGInsertions(&data)
		sources := GetSources(ANIMAL | BACTERIA)
		fmt.Println("CGGCGG")
		CountPattern(sources, []byte("CGGCGG"))

		fmt.Println("FCS")
		CountPattern(sources, []byte("CTCCTCGGCGGG"))
	}

	if countFCSAlternatives {
		sources := GetSources(BACTERIA)
		for i, _ := range sources {
			source := &sources[i]
			CountFCSAlternatives(nil, source)
		}
	}

	if fcsAlternatives {
		CheckAlternatives()
	}

	if outputGob {
		data.OutputMatches("matches.txt")
		data.OutputInsertions("insertion-data.txt")
	}
}
