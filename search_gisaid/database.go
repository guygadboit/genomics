package main

import (
	"bufio"
	"encoding/gob"
	"flag"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
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

type Range struct {
	Start, End int
}

func ParseRanges(s string) []Range {
	ret := make([]Range, 0)
	if s == "" {
		return ret
	}
	fields := strings.Split(s, ",")
	for _, f := range fields {
		subFields := strings.Split(f, "-")

		start := utils.Atoi(subFields[0])
		var end int

		switch len(subFields) {
		case 1:
			end = start
		case 2:
			end = utils.Atoi(subFields[1])
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
	GisaidAccession   string       // 0
	Isolate           string       // 1
	SubmissionDate    Date         // 2
	CollectionDate    Date         // 3
	PangolinLineage   string       // 4
	Country           string       // 5
	Region            string       // 6
	City              string       // 7
	Length            int          // 8
	Host              string       // 9
	Divergence        int          // 10
	NucleotideChanges []Mutation   // 11
	Deletions         []Range      // 12
	Insertions        []Insertion  // 13
	AAChanges         []AAMutation // 14
	WhoClade          string       // 15
	NextstrainClade   string       // 16
	Continent         string       // 17
	ToBeExcluded      int          // 18
}

func (r *Record) ToString() string {
	return fmt.Sprintf("%s %s %s %s %s", r.GisaidAccession,
		r.CollectionDate.ToString(), r.Country, r.Region, r.City)
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
	r.Deletions = ParseRanges(fields[12])
	r.Insertions = ParseInsertions(fields[13])
	r.AAChanges = ParseAAMutations(fields[14])
	r.WhoClade = fields[15]
	r.NextstrainClade = fields[16]
	r.Continent = fields[17]
	r.ToBeExcluded = Atoi(fields[18])
}

type Database struct {
	Records       []Record
	MutationIndex map[int][]int
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

func SantaCatarina(db Database) {
	muts := ParseMutations("A5706G,C14408T")
	matches := db.Search(muts)

	for k, _ := range matches {
		r := db.Records[k]
		fmt.Printf("%s %s %s %s\n",
			r.Country, r.Region, r.City, r.CollectionDate.ToString())
	}

	for _, r := range db.Records {
		if r.Region == "Santa Catarina" {
			fmt.Println(r.ToString())
		}
	}
}

func main() {
	var save bool

	flag.BoolVar(&save, "s", false, "Parse and save to gob")
	flag.Parse()

	var database Database

	if save {
		database = Parse("/fs/h/genomes/GISAID/gisaid2020.tsv")
		fmt.Printf("Parsed %d records\n", len(database.Records))
		database.BuildMutationIndex()
		fmt.Printf("Built index\n")
		database.Save("db.gob")
		return
	}

	database.Load("db.gob")
	fmt.Printf("Loaded\n")

	// Get rid of this later if you want to make this more generic
	SantaCatarina(database)
}
