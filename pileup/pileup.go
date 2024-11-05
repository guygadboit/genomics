package pileup

import (
	"genomics/utils"
	"slices"
	"strings"
)

type Read struct {
	Nt    byte
	Depth int
}

type Record struct {
	Pos   int
	Reads []Read // Sorted by highest depth first
}

type Pileup struct {
	// The index maps genome positions to positions in the reads array. It will
	// just be 1:1 most of the time since you will have reads at every position
	Index   map[int]int
	Records []Record
}

func (p *Pileup) Init() {
	p.Index = make(map[int]int)
	p.Records = make([]Record, 0)
}

func (p *Pileup) Add(pos int, reads []Read) {
	p.Records = append(p.Records, Record{pos, reads})
	p.Index[pos] = len(p.Records)-1
}

func (p *Pileup) Get(pos int) Record {
	return p.Records[p.Index[pos]]
}

type counter map[byte]int

func parseReadBases(s string) []Read {
	s = strings.ToUpper(s)
	counts := make(counter)
	for _, c := range []byte(s) {
		switch c {
		case 'A':
			fallthrough
		case 'G':
			fallthrough
		case 'T':
			fallthrough
		case 'C':
			counts[c]++
		}
	}

	ret := make([]Read, 0, len(counts))
	for k, v := range counts {
		ret = append(ret, Read{k, v})
	}

	slices.SortFunc(ret, func(a, b Read) int {
		if a.Depth < b.Depth {
			return 1
		}
		if a.Depth > b.Depth {
			return -1
		}
		return 0
	})

	return ret
}

func Parse(fname string) (Pileup, error) {
	var err error
	var ret Pileup
	ret.Init()

	utils.Lines(fname, func(line string, lineErr error) bool {
		if lineErr != nil {
			err = lineErr
			return false
		}
		fields := strings.Split(line, "\t")

		pos := utils.Atoi(fields[1]) - 1
		reads := parseReadBases(fields[4])

		if len(reads) != 0 {
			ret.Add(pos, reads)
		}
		return true
	})

	if err != nil {
		return ret, err
	}

	return ret, nil
}
