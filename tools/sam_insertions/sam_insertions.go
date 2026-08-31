package main

import (
	"flag"
	"fmt"
	"genomics/utils"
	"regexp"
	"strings"
)

type Insertion struct {
	Pos utils.OneBasedPos
	Nts []byte
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

func FindInsertions(samLine string) []Insertion {
	ret := make([]Insertion, 0)
	fields := strings.Fields(samLine)
	if len(fields) < 10 {
		return ret
	}

	nts := fields[9]
	cigar := ParseCigar(fields[5])

	pos := 0
	for _, cf := range cigar {
		switch cf.Op {
		case 'M':
			fallthrough
		case 'S':
			pos += cf.Value
		case 'I':

			// FIXME YOU ARE HERE

			
		}

	}


	fmt.Println(cigar, cf, nts)
	return ret
}

func main() {
	flag.Parse()
	fname := flag.Arg(0)

	utils.Lines(fname, func(line string, err error) bool {
		FindInsertions(line)
		return true
	})
}
