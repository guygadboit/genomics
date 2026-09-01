package main

import (
	"flag"
	"fmt"
	"genomics/utils"
	"log"
	"path"
	"regexp"
	"strings"
	"os"
	"io"
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

func FindInsertions(samLine string, minLen, maxLen int) []Insertion {
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
			if maxLen != -1 && cf.Value > maxLen {
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

type CountedInsertion struct {
	Insertion
	Count int
}

func ShowCountedInsertions(insMap map[string]CountedInsertion, fp io.Writer) {
	inss := make([]CountedInsertion, 0, len(insMap))
	for _, ci := range insMap {
		inss = append(inss, ci)
	}
	utils.SortByKey(inss, false, func(ci CountedInsertion) int {
		return len(ci.Insertion.Nts)
	})
	for i, ins := range inss {
		fmt.Fprintf(fp, "%d ins_%d:%s (%d reads)\n",
			i+1, ins.Pos, string(ins.Nts), ins.Count)
	}
}

func main() {
	var (
		minLen, maxLen int
		gzip bool
	)

	flag.IntVar(&minLen, "min", 12, "Minimum length")
	flag.IntVar(&maxLen, "max", -1, "Maximum length")
	flag.BoolVar(&gzip, "z", false, "Gzip output")

	flag.Parse()

	var writeFn func(string) (*os.File, io.WriteCloser)
	if gzip {
		writeFn = utils.WriteFileGz
	} else {
		writeFn = utils.WriteFile
	}

	for _, fname := range flag.Args() {

		allInsertions := make(map[string]CountedInsertion)

		var lineNo int
		utils.Lines(fname, func(line string, err error) bool {
			lineNo++
			if line[0] == '@' {
				return true // skip header lines
			}
			insertions := FindInsertions(line, minLen, maxLen)
			for _, ins := range insertions {
				key := ins.String()
				record, there := allInsertions[key]
				if there {
					allInsertions[key] = CountedInsertion{ins, record.Count + 1}
				} else {
					allInsertions[key] = CountedInsertion{ins, 1}
				}
			}
			return true
		})

		_, f := path.Split(fname)
		outFname := utils.BaseName(f) + ".ins"
		if gzip {
			outFname += ".gz"
		}

		fd, fp := writeFn(outFname)
		defer fd.Close()

		ShowCountedInsertions(allInsertions, fp)
		fp.Close()
		fmt.Printf("Wrote %s\n", outFname)
	}
}
