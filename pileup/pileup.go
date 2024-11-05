package pileup

import (
	"genomics/utils"
	"strings"
	"fmt"
)

type Record struct {
	Pos   int
	Nt    byte
	Depth int
}

type Pileup []Record

type Counter map[byte]int

func (c Counter) findMax() (byte, int, int) {
	var bestByte byte
	var bestCount int

	for k, v := range c {
		if v > bestCount {
			bestByte = k
			bestCount = v
		}
	}

	var numTie int
	for _, v := range c {
		if v == bestCount {
			numTie++
		}
	}

	return bestByte, bestCount, numTie
}

func parseReadBases(s string) (byte, int, int) {
	s = strings.ToUpper(s)
	counts := make(Counter)
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
	return counts.findMax()
}

func Parse(fname string) (Pileup, error) {
	var err error
	ret := make(Pileup, 0)

	utils.Lines(fname, func(line string, lineErr error) bool {
		if lineErr != nil {
			err = lineErr
			return false
		}
		fields := strings.Split(line, "\t")

		pos := utils.Atoi(fields[1]) - 1
		nt, depth, numTie := parseReadBases(fields[4])

		if nt != 0 {
			record := Record{pos, nt, depth}
			ret = append(ret, record)
			if numTie > 1 {
				fmt.Printf("%d-way ambiguity at %d\n", numTie, pos+1)
			}
		}
		return true
	})

	if err != nil {
		return nil, err
	}

	return ret, nil
}
