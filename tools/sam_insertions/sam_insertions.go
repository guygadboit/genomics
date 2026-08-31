package main

import (
	"flag"
	"fmt"
	"genomics/utils"
	"log"
	"regexp"
	"strings"
)

type Insertion struct {
	Pos utils.OneBasedPos
	Nts []byte
}

func (ins Insertion) String() string {
	return fmt.Sprintf("%d:%s", ins.Pos, string(ins.Nts))
}

type CigarField struct {
	Op    byte
	Value int
}

func (cf CigarField) String() string {
	return fmt.Sprintf("%d%c", cf.Value, cf.Op)
}

type Cigar []CigarField

func ParseCigar(cigar string) Cigar {
	pat := regexp.MustCompile(`\d+[A-Z]`)
	matches := pat.FindAllString(cigar, -1)
	ret := make(Cigar, len(matches))
	for i, m := range matches {
		n := len(m)
		ret[i] = CigarField{m[n-1], utils.Atoi(m[:n-1])}
	}
	return ret
}

func FindInsertions(samLine string, minLen int) []Insertion {
	ret := make([]Insertion, 0)
	fields := strings.Fields(samLine)
	if len(fields) < 10 {
		return ret
	}

	pos := utils.OneBasedPos(utils.Atoi(fields[3]))
	nts := []byte(fields[9])
	cigar := ParseCigar(fields[5])

	readPos := 0
	for _, cf := range cigar {
		switch cf.Op {
		case 'S':
			fallthrough
		case 'M':
			readPos += cf.Value
		case 'N':
			fallthrough
		case 'D':
			pos += utils.OneBasedPos(cf.Value)
		case 'I':
			if cf.Value < minLen {
				continue
			}
			end := readPos + cf.Value
			if end < len(nts) {
				// -1 here seems to be the convention for how insertions are
				// written down (verified with seqkit mutate -i)
				ins := Insertion{pos + utils.OneBasedPos(readPos) - 1,
					nts[readPos:end]}
				ret = append(ret, ins)
				// fmt.Printf("%d: %s\n", ins.Pos, string(ins.Nts))
			} else {
				// This shouldn't happen unless your SAM file is invalid
				log.Printf("Out of bounds!")
			}
		}
	}

	return ret
}

func main() {
	flag.Parse()
	fname := flag.Arg(0)

	type countedInsertion struct {
		Insertion
		Count int
	}
	allInsertions := make(map[string]countedInsertion)

	var lineNo int
	utils.Lines(fname, func(line string, err error) bool {
		lineNo++
		if line[0] == '@' {
			return true // skip header lines
		}
		insertions := FindInsertions(line, 12)
		for _, ins := range insertions {
			key := ins.String()
			record, there := allInsertions[key]
			if there {
				allInsertions[key] = countedInsertion{ins, record.Count + 1}
			} else {
				allInsertions[key] = countedInsertion{ins, 1}
			}
		}
		return true
	})
	for k, v := range allInsertions {
		fmt.Printf("%s: %d\n", k, v.Count)
	}
}
