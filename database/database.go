package database

import (
	"bufio"
	"encoding/gob"
	"errors"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"io"
	"log"
	"os"
	"reflect"
	"slices"
	"strings"
	"time"
)

const (
	ROOT     = "/fs/f/genomes/GISAID/"
	GOB_NAME = ROOT + "gisaid2020.gob"
)

// These are actually just indexes into Database.Records
type Id int

type Silence int

const (
	UNKNOWN Silence = iota
	SILENT
	NON_SILENT
	NOT_IN_ORF
)

type Mutation struct {
	Pos     utils.OneBasedPos
	From    byte
	To      byte
	Silence Silence
}

type Mutations []Mutation
type AAMutations []AAMutation

func (m *Mutation) ToString() string {
	var silence string
	switch m.Silence {
	case SILENT:
		silence = "*"
	case NOT_IN_ORF:
		silence = "@"
	}
	return fmt.Sprintf("%c%d%c%s", m.From, m.Pos, m.To, silence)
}

func (muts Mutations) ToString() string {
	s := make([]string, len(muts))
	for i, m := range muts {
		s[i] = m.ToString()
	}
	return strings.Join(s, ",")
}

type AAMutation struct {
	Mutation
	Gene string
}

func (m *AAMutation) ToString() string {
	return fmt.Sprintf("%s:%c%d%c", m.Gene, m.From, m.Pos, m.To)
}

func (muts AAMutations) ToString() string {
	s := make([]string, len(muts))
	for i, m := range muts {
		s[i] = m.ToString()
	}
	return strings.Join(s, ",")
}

// Parses a string like "A17858G,C17747T". No spaces allowed!
func ParseMutations(s string) []Mutation {
	ret := make([]Mutation, 0)
	if s == "" {
		return ret
	}
	fields := strings.Split(s, ",")
	for _, f := range fields {
		n := len(f)
		mut := Mutation{utils.OneBasedPos(utils.Atoi(f[1 : n-1])),
			f[0], f[n-1], UNKNOWN}
		ret = append(ret, mut)
	}
	return ret
}

func ParseAAMutations(s string) []AAMutation {
	ret := make([]AAMutation, 0)
	if s == "" {
		return ret
	}

	fields := strings.Split(s, ",")
	for _, field := range fields {
		subFields := strings.Split(field, ":")
		gene := subFields[0]
		f := subFields[1]
		n := len(f)
		mut := AAMutation{Mutation{utils.OneBasedPos(utils.Atoi(f[1 : n-1])),
			f[0], f[n-1], UNKNOWN}, gene}
		ret = append(ret, mut)
	}
	return ret
}

type Range struct {
	Start, End utils.OneBasedPos
}

func (r *Range) ToString() string {
	return fmt.Sprintf("%d-%d", r.Start, r.End)
}

func ParseRanges(s string) []Range {
	ret := make([]Range, 0)
	if s == "" {
		return ret
	}
	fields := strings.Split(s, ",")
	for _, f := range fields {
		subFields := strings.Split(f, "-")

		start := utils.OneBasedPos(utils.Atoi(subFields[0]))
		var end utils.OneBasedPos

		switch len(subFields) {
		case 1:
			end = start
		case 2:
			end = utils.OneBasedPos(utils.Atoi(subFields[1]))
		default:
			log.Fatalf("Invalid range %s\n", f)
		}
		ret = append(ret, Range{start, end})
	}
	return ret
}

type Insertion struct {
	Pos      int
	Sequence []byte
}

func ParseInsertions(s string) []Insertion {
	ret := make([]Insertion, 0)
	if s == "" {
		return ret
	}
	fields := strings.Split(s, ",")
	for _, f := range fields {
		subFields := strings.Split(f, ":")
		pos := utils.Atoi(subFields[0])
		seq := []byte(subFields[1])
		ret = append(ret, Insertion{pos, seq})
	}
	return ret
}

type Record struct {
	Id                Id
	GisaidAccession   string      // 0
	Isolate           string      // 1
	SubmissionDate    time.Time   // 2
	CollectionDate    time.Time   // 3
	PangolinLineage   string      // 4
	Country           string      // 5
	Region            string      // 6
	City              string      // 7
	Length            int         // 8
	Host              string      // 9
	Divergence        int         // 10
	NucleotideChanges Mutations   // 11
	Deletions         []Range     // 12
	Insertions        []Insertion // 13
	AAChanges         AAMutations // 14
	WhoClade          string      // 15
	NextstrainClade   string      // 16
	Continent         string      // 17
	ToBeExcluded      int         // 18
}

func (r *Record) ToString() string {
	return fmt.Sprintf("%s %s %s %s %s", r.GisaidAccession,
		r.CollectionDate.Format(time.DateOnly), r.Country, r.Region, r.City)
}

func (r *Record) FilterNucleotideChanges(silence ...Silence) Mutations {
	ret := make(Mutations, 0)
	for _, c := range r.NucleotideChanges {
		for _, s := range silence {
			if c.Silence == s {
				ret = append(ret, c)
			}
		}
	}
	return ret
}

// A more detailed summary
func (r *Record) Summary() string {
	return fmt.Sprintf("%s %s %s %s %s: %s %s; %d %d",
		r.GisaidAccession,
		r.CollectionDate.Format(time.DateOnly),
		r.Country,
		r.Region,
		r.City,
		r.NucleotideChanges.ToString(),
		r.AAChanges.ToString(),
		len(r.Insertions),
		len(r.Deletions))
}

func (r *Record) DeletionsSummary() string {
	s := make([]string, len(r.Deletions))
	for i, d := range r.Deletions {
		s[i] = fmt.Sprintf("D%s", d.ToString())
	}
	return strings.Join(s, ",")
}

func Atoi(s string) int {
	if s == "" {
		return 0
	}
	return utils.Atoi(s)
}

func (r *Record) Parse(line string) {
	fields := strings.Split(line, "\t")
	r.GisaidAccession = fields[0]
	r.Isolate = fields[1]
	r.SubmissionDate, _ = time.Parse(time.DateOnly, fields[2])
	r.CollectionDate, _ = time.Parse(time.DateOnly, fields[3])
	r.PangolinLineage = fields[4]
	r.Country = fields[5]
	r.Region = fields[6]
	r.City = fields[7]
	r.Length = Atoi(fields[8])
	r.Host = fields[9]
	r.Divergence = Atoi(fields[10])
	r.NucleotideChanges = ParseMutations(fields[11])
	r.Deletions = ParseRanges(fields[12])
	r.Insertions = ParseInsertions(fields[13])
	r.AAChanges = ParseAAMutations(fields[14])
	r.WhoClade = fields[15]
	r.NextstrainClade = fields[16]
	r.Continent = fields[17]
	r.ToBeExcluded = Atoi(fields[18])
}

type Database struct {
	Records        []Record
	MutationIndex  map[utils.OneBasedPos]IdSet
	AccessionIndex map[string]Id
}

func (d *Database) Init() {
	d.Records = make([]Record, 0)
}

func (d *Database) Add(r *Record) {
	r.Id = Id(len(d.Records))
	d.Records = append(d.Records, *r)
}

func (d *Database) Get(id Id) *Record {
	return &d.Records[id]
}

func (d *Database) BuildMutationIndex() {
	d.MutationIndex = make(map[utils.OneBasedPos]IdSet)

	for i, r := range d.Records {
		for _, mut := range r.NucleotideChanges {
			if d.MutationIndex[mut.Pos] == nil {
				d.MutationIndex[mut.Pos] = make(IdSet)
			}
			d.MutationIndex[mut.Pos][Id(i)] = true
		}
	}
}

func (d *Database) BuildAccessionIndex() {
	d.AccessionIndex = make(map[string]Id)

	for i, r := range d.Records {
		d.AccessionIndex[r.GisaidAccession] = Id(i)
	}
}

type IdSet map[Id]bool

// minMatches is the minimum number of matches-- so 1 for "Any" or len(muts)
// for "All".
func (d *Database) SearchByMutPosition(pos []utils.OneBasedPos,
	minMatches int) IdSet {
	matches := make(map[Id]int) // How many matches for each record?

	for _, pos := range pos {
		found, there := d.MutationIndex[pos]
		if there {
			for k, _ := range found {
				matches[k]++
			}
		}
	}

	if minMatches > len(pos) {
		minMatches = len(pos)
	}

	ret := make(IdSet)
	for k, v := range matches {
		if v >= minMatches {
			ret[k] = true
		}
	}

	return ret
}

func (d *Database) GetByAccession(accNum ...string) []Id {
	ret := make([]Id, 0)

	for _, an := range accNum {
		id, there := d.AccessionIndex[an]
		if there {
			ret = append(ret, id)
		}
	}
	return ret
}

func (d *Database) Save(fname string) {
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	enc := gob.NewEncoder(fp)
	err = enc.Encode(d)
	if err != nil {
		log.Fatal(err)
	}
	fp.Flush()
}

func (d *Database) Load(fname string) {
	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)
	dec := gob.NewDecoder(fp)
	err = dec.Decode(d)
	if err != nil {
		log.Fatal(err)
	}
}

func (d *Database) Parse(fname string) {
	d.Init()
	fp := utils.NewFileReader(fname)
	defer fp.Close()

loop:
	for count := 0; ; count++ {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		if count == 0 { // skip the headings line
			continue
		}

		line = strings.TrimSpace(line)

		var record Record
		record.Parse(line)
		d.Add(&record)
	}
	d.BuildMutationIndex()
	d.BuildAccessionIndex()
}

// Return whichever muts in muts this record has
func (r *Record) HasMuts(muts Mutations) Mutations {
	ret := make([]Mutation, 0)
	for _, s := range muts {
		for _, m := range r.NucleotideChanges {
			if reflect.DeepEqual(s, m) {
				ret = append(ret, m)
			}
		}
	}
	return ret
}

// You could use generics for this
func (r *Record) HasAAMuts(muts AAMutations) AAMutations {
	ret := make([]AAMutation, 0)
	for _, s := range muts {
		for _, m := range r.AAChanges {
			if reflect.DeepEqual(s, m) {
				ret = append(ret, m)
			}
		}
	}
	return ret
}

type Key int

const (
	COLLECTION_DATE Key = iota
	COUNTRY
	REGION
)

func (d *Database) Sort(ids []Id, key Key) {
	var cmp func(a, b *Record) int

	switch key {
	case COLLECTION_DATE:
		cmp = func(a, b *Record) int {
			return a.CollectionDate.Compare(b.CollectionDate)
		}
	default:
		log.Fatal("Sort key not implemented")
	}

	slices.SortFunc(ids, func(a, b Id) int {
		return cmp(&d.Records[a], &d.Records[b])
	})
}

// nil for ids means everything
func (d *Database) Filter(ids IdSet, fun func(r *Record) bool) IdSet {
	ret := make(IdSet)

	if ids == nil {
		for _, r := range d.Records {
			if fun(&r) {
				ret[r.Id] = true
			}
		}
	} else {
		for id, _ := range ids {
			if fun(&d.Records[id]) {
				ret[id] = true
			}
		}
	}
	return ret
}

func (d *Database) DetermineSilence(reference *genomes.Genomes) {
	cache := make(map[Mutation]Silence)

	determine := func(mut Mutation) Silence {
		silence, there := cache[mut]
		if !there {
			isSilent, _, err := genomes.IsSilentWithReplacement(reference,
				int(mut.Pos)-1, 0, 0, []byte{mut.To})
			if err != nil {
				silence = NOT_IN_ORF
			} else if isSilent {
				silence = SILENT
			} else {
				silence = NON_SILENT
			}
			cache[mut] = silence
		}
		return silence
	}

	for i, r := range d.Records {
		for j, mut := range r.NucleotideChanges {
			d.Records[i].NucleotideChanges[j].Silence = determine(mut)
		}
	}
}

/*
Use the first genome in reference, and output it and the second one in an
alignment.
*/
func (d *Database) Reconstruct(id Id,
	reference *genomes.Genomes, name string) (*genomes.Genomes, error) {

	ret := reference.Filter(0, 0)
	ret.DeepCopy(1)
	record := &d.Records[id]

	for _, mut := range record.NucleotideChanges {
		pos := mut.Pos - 1
		if ret.Nts[1][pos] != mut.From {
			return nil, errors.New("Wrong reference")
		}

		ret.Nts[1][pos] = mut.To
	}

	for _, ins := range record.Insertions {
		pos := ins.Pos - 1
		n := len(ins.Sequence)
		fmt.Printf("Insertion at %d length %d\n", pos, n)

		ret.Nts[0] = utils.Insert(ret.Nts[0], pos, n)
		ret.Nts[1] = utils.Insert(ret.Nts[1], pos, n)

		copy(ret.Nts[1][pos:pos+n], ins.Sequence)
		for i := pos; i < pos+n; i++ {
			ret.Nts[0][i] = '-'
		}
	}

	for _, del := range record.Deletions {
		fmt.Printf("Deletion at %d-%d\n", del.Start-1, del.End)
		for i := del.Start - 1; i < del.End; i++ {
			ret.Nts[1][i] = '-'
		}
	}

	ret.Names[1] = name
	return ret, nil
}

func NewDatabase() *Database {
	var ret Database
	ret.Load(GOB_NAME)
	return &ret
}
