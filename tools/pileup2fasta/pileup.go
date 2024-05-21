package main

import (
	"genomics/utils"
	"strings"
)

type Record struct {
	pos   int
	nt    byte
	depth int
}

type Pileup []Record

type Counter map[byte]int

func (c Counter) findMax() (byte, int) {
	var bestByte byte
	var bestCount int

	for k, v := range c {
		if v > bestCount {
			bestByte = k
			bestCount = v
		}
	}
	return bestByte, bestCount
}

func parseReadBases(s string) (byte, int) {
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
		nt, depth := parseReadBases(fields[4])

		if nt != 0 {
			record := Record{pos, nt, depth}
			ret = append(ret, record)
		}
		return true
	})

	if err != nil {
		return nil, err
	}

	return ret, nil
}
