package main

import (
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"strings"
)

type Date struct {
	Y, M, D int
}

func (d *Date) Parse(s string) {
	fields := strings.Split(s, "-")
	d.Y = utils.Atoi(fields[0])

	if len(fields) >= 2 {
		d.M = utils.Atoi(fields[1])
	}
	if len(fields) == 3 {
		d.D = utils.Atoi(fields[2])
	}
}

func (d *Date) ToString() string {
	return fmt.Sprintf("%d-%02d-%02d", d.Y, d.M, d.D)
}

type Mutation struct {
	Pos  int
	From byte
	To   byte
}

type AAMutation struct {
	Mutation
	Gene string
}

func ParseMutations(s string) []Mutation {
	ret := make([]Mutation, 0)
	if s == "" {
		return ret
	}
	fields := strings.Split(s, ",")
	for _, f := range fields {
		n := len(f)
		mut := Mutation{utils.Atoi(f[1 : n-1]), f[0], f[n-1]}
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
		mut := AAMutation{Mutation{utils.Atoi(f[1 : n-1]),
			f[0], f[n-1]}, gene}
		ret = append(ret, mut)
	}
	return ret
}

func ParseInts(s string) []int {
	ret := make([]int, 0)
	fields := strings.Split(s, ",")
	for _, f := range fields {
		ret = append(ret, utils.Atoi(f))
	}
	return ret
}

type Record struct {
	GisaidAccession   string     // 0
	Isolate           string     // 1
	SubmissionDate    Date       // 2
	CollectionDate    Date       // 3
	PangolinLineage   string     // 4
	Country           string     // 5
	Region            string     // 6
	City              string     // 7
	Length            int        // 8
	Host              string     // 9
	Divergence        int        // 10
	NucleotideChanges []Mutation // 11

	// These seem to be ranges or something, not just lists. So let's just keep
	// them as strings for now.
	// Deletions         []int        // 12
	// Insertions        []int        // 13
	Deletions  string // 12
	Insertions string // 13

	AAChanges       []AAMutation // 14
	WhoClade        string       // 15
	NextstrainClade string       // 16
	Continent       string       // 17
	ToBeExcluded    int          // 18
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
	r.SubmissionDate.Parse(fields[2])
	r.CollectionDate.Parse(fields[3])
	r.PangolinLineage = fields[4]
	r.Country = fields[5]
	r.Region = fields[6]
	r.City = fields[7]
	r.Length = Atoi(fields[8])
	r.Host = fields[9]
	r.Divergence = Atoi(fields[10])
	r.NucleotideChanges = ParseMutations(fields[11])
	// r.Deletions = ParseInts(fields[12])
	r.Deletions = fields[12]
	// r.Insertions = ParseInts(fields[13])
	r.Insertions = fields[13]
	r.AAChanges = ParseAAMutations(fields[14])
	r.WhoClade = fields[15]
	r.NextstrainClade = fields[16]
	r.Continent = fields[17]
	r.ToBeExcluded = Atoi(fields[18])
}

type Database struct {
	Records	[]Record
	MutationIndex	map[int][]int
}

func (d *Database) Init() {
	d.Records = make([]Record, 0)
}

func (d *Database) Add(r *Record) {
	d.Records = append(d.Records, *r)
}

func (d *Database) BuildMutationIndex() {
	d.MutationIndex = make(map[int][]int)

	for i, r := range d.Records {
		for _, mut := range r.NucleotideChanges {
			d.MutationIndex[mut.Pos] = append(d.MutationIndex[mut.Pos], i)
		}
	}
}

type Set map[int]bool

func Intersection(a, b Set) Set {
	ret := make(Set)
	for k, _ := range a {
		if b[k] {
			ret[k] = true
		}
	}
	return ret
}

func (d *Database) Search(muts []Mutation) Set {
	var candidates Set
	for _, mut := range muts {
		found, there := d.MutationIndex[mut.Pos]
		if there {
			s := utils.ToSet(found)
			if candidates == nil {
				candidates = s
			} else {
				candidates = Intersection(candidates, s)
			}
		}
	}
	return candidates
}

func Parse(fname string) Database {
	var ret Database
	ret.Init()
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
		ret.Add(&record)
	}
	return ret
}

func main() {
	database := Parse("/fs/h/genomes/GISAID/gisaid2020.tsv")
	fmt.Printf("Parsed %d records\n", len(database.Records))
	database.BuildMutationIndex()
	fmt.Printf("Built index\n")
	// FIXME you could gob it now

	muts := ParseMutations("A5706G,C14408T")
	matches := database.Search(muts)

	for k, _ := range matches {
		r := database.Records[k]
		fmt.Printf("%s %s %s %s\n",
			r.Country, r.Region, r.City, r.CollectionDate.ToString())
	}
}
