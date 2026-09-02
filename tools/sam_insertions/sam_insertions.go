package main

import (
	"flag"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"path"
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

func (im InsertionMap) Add(ins *Insertion, count int) {
	key := ins.String()
	record, there := im[key]
	if there {
		im[key] = CountedInsertion{*ins, record.Count + count}
	} else {
		im[key] = CountedInsertion{*ins, count}
	}
}

type InsertionMap map[string]CountedInsertion

func (im InsertionMap) Combine(other InsertionMap) {
	for _, v := range other {
		im.Add(&v.Insertion, 1)
	}
}

func (insMap InsertionMap) Display(fp io.Writer) {
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

func (insMap InsertionMap) FromFile(fname string) {
	pat := regexp.MustCompile(`\d+ ins_(\d+):([GNATC]+) \((\d+) reads\)`)
	utils.Lines(fname, func(line string, err error) bool {
		groups := pat.FindAllStringSubmatch(line, -1)
		ins := Insertion{
			utils.OneBasedPos(utils.Atoi(groups[0][1])),
			[]byte(groups[0][2])}
		count := utils.Atoi(groups[0][3])
		insMap.Add(&ins, count)
		return true
	})
}

func Load(fname string, minLen, maxLen int) InsertionMap {
	ret := make(InsertionMap)

	var lineNo int
	utils.Lines(fname, func(line string, err error) bool {
		lineNo++
		if line[0] == '@' {
			return true // skip header lines
		}
		insertions := FindInsertions(line, minLen, maxLen)
		for _, ins := range insertions {
			key := ins.String()
			record, there := ret[key]
			if there {
				ret[key] = CountedInsertion{ins, record.Count + 1}
			} else {
				ret[key] = CountedInsertion{ins, 1}
			}
		}
		return true
	})
	return ret
}

func Merge(outName string, fnames ...string) {
	all := make(InsertionMap)

	for _, fname := range fnames {
		im := make(InsertionMap)
		im.FromFile(fname)
		all.Combine(im)
	}

	fd, fp := utils.WriteFileGz(outName)
	defer fd.Close()

	all.Display(fp)
	fp.Close()
}

func main() {
	var (
		minLen, maxLen int
		gzip           bool
		outName        string
		merge          bool
	)

	flag.IntVar(&minLen, "min", 12, "Minimum length")
	flag.IntVar(&maxLen, "max", -1, "Maximum length")
	flag.BoolVar(&gzip, "z", false, "Gzip output")
	flag.StringVar(&outName, "0", "output.ins.gz", "Output name for merge")
	flag.BoolVar(&merge, "merge", false, "Merge a bunch of .ins.gz files")

	flag.Parse()

	if merge {
		Merge(outName, flag.Args()...)
		return
	}

	var writeFn func(string) (*os.File, io.WriteCloser)
	if gzip {
		writeFn = utils.WriteFileGz
	} else {
		writeFn = utils.WriteFile
	}

	for _, fname := range flag.Args() {
		allInsertions := Load(fname, minLen, maxLen)

		_, f := path.Split(fname)
		outFname := utils.BaseName(f) + ".ins"
		if gzip {
			outFname += ".gz"
		}

		fd, fp := writeFn(outFname)
		defer fd.Close()

		allInsertions.Display(fp)
		fp.Close()
		fmt.Printf("Wrote %s\n", outFname)
	}
}
