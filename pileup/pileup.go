package pileup

import (
	"fmt"
	"genomics/utils"
	"slices"
	"strings"
)

type Read struct {
	Nt    byte
	Depth int
}

type Record struct {
	Pos        int
	Reads      []Read // Sorted by highest depth first
	TotalDepth int    // Summed over all the reads
}

func (r *Record) GetDepthOf(nt byte) int {
	if r == nil {
		return 0
	}
	for _, read := range r.Reads {
		if read.Nt == nt {
			return read.Depth
		}
	}
	return 0
}

type Pileup struct {
	// The index maps genome positions to positions in the reads array. It will
	// just be 1:1 most of the time since you will have reads at every position
	Index   map[int]int
	Records []Record
	MaxPos  int
}

func (p *Pileup) Init() {
	p.Index = make(map[int]int)
	p.Records = make([]Record, 0)
}

// reads are expected to already be sorted by highest depth first
func (p *Pileup) Add(pos int, reads []Read) {
	var totalDepth int
	for _, read := range reads {
		totalDepth += read.Depth
	}
	p.Records = append(p.Records, Record{pos, reads, totalDepth})
	p.Index[pos] = len(p.Records) - 1
	if pos > p.MaxPos {
		p.MaxPos = pos
	}
}

func CombineReads(a, b []Read) []Read {
	m := make(map[byte]int)
	for _, r := range a {
		m[r.Nt] += r.Depth
	}
	for _, r := range b {
		m[r.Nt] += r.Depth
	}
	ret := make([]Read, 0, len(m))
	for k, v := range m {
		ret = append(ret, Read{k, v})
	}
	utils.SortByKey(ret, false, func(r Read) int { return r.Depth })
	return ret
}

func CloneReads(r []Read) []Read {
	ret := make([]Read, len(r))
	copy(ret, r)
	return ret
}

// Return a new Pileup which is the sum of two others
func Combine(a, b *Pileup) *Pileup {
	var ret Pileup
	ret.Init()

	for pos := 0; pos < max(a.MaxPos, b.MaxPos); pos++ {
		aRec := a.Get(pos)
		bRec := b.Get(pos)

		var reads []Read
		if aRec != nil && bRec != nil {
			reads = CombineReads(aRec.Reads, bRec.Reads)
		} else if aRec != nil {
			reads = CloneReads(aRec.Reads)
		} else if bRec != nil {
			reads = CloneReads(bRec.Reads)
		}

		ret.Add(pos, reads)
	}
	return &ret
}

func (p *Pileup) Get(pos int) *Record {
	recPos, there := p.Index[pos]
	if !there {
		return nil
	}
	return &p.Records[recPos]
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

/*
Parse the output of samtools mpileup
*/
func Parse(fname string) (*Pileup, error) {
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
		return nil, err
	}

	return &ret, nil
}

func FormatRecord(record *Record) string {
	items := make([]string, len(record.Reads))
	for i, read := range record.Reads {
		items[i] = fmt.Sprintf("%cx%d", read.Nt, read.Depth)
	}
	return fmt.Sprintf("%d: %s", record.Pos+1, strings.Join(items, ", "))
}

func (pu *Pileup) Show(onlyPos []int) {
	displayRecord := func(pos int) {
		recordI, there := pu.Index[pos]
		if !there {
			fmt.Printf("%d:\n", pos+1)
			return
		}
		record := &pu.Records[recordI]
		fmt.Println(FormatRecord(record))
	}

	if onlyPos != nil {
		for _, pos := range onlyPos {
			displayRecord(pos)
		}
	} else {
		for i := 0; i <= pu.MaxPos; i++ {
			displayRecord(i)
		}
	}
}

func (p *Pileup) Coverage(minDepth int) float64 {
	var got, total float64
	for i := 0; i <= p.MaxPos; i++ {
		r := p.Get(i)
		total++
		if r == nil {
			continue
		}
		if r.Reads[0].Depth < minDepth {
			continue
		}
		got++
	}
	return got / total
}

/*
Parse the output of our own "show" option back in
*/
func Parse2(fname string) (*Pileup, error) {
	var ret Pileup
	ret.Init()

	utils.Lines(fname, func(line string, lineErr error) bool {
		line = strings.TrimRight(line, " ")
		fields := strings.Split(line, ":")
		pos := utils.Atoi(fields[0]) - 1

		if fields[1] == "" {
			return true
		}

		reads := make([]Read, 0)
		readData := strings.Split(fields[1], ",")
		for _, d := range readData {
			fields := strings.Split(d, "x")
			nt, depth := strings.Trim(fields[0], " "), utils.Atoi(fields[1])
			reads = append(reads, Read{nt[0], depth})
		}
		ret.Add(pos, reads)
		return true
	})
	return &ret, nil
}
