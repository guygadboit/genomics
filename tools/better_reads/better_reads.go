package main

import (
	"flag"
	"fmt"
	"genomics/utils"
	"strings"
)

type Read struct {
	Name  string
	Score int
}

type ReadSet struct {
	Name      string
	Reads     []Read
	NameIndex map[string]int
}

func (r *ReadSet) BuildNameIndex() {
	r.NameIndex = make(map[string]int)
	for i, read := range r.Reads {
		r.NameIndex[read.Name] = i
	}
}

func (r *ReadSet) GetByName(name string) *Read {
	i, there := r.NameIndex[name]
	if !there {
		return nil
	}
	return &r.Reads[i]
}

func (r *ReadSet) Better(others *ReadSet, verbose bool) *ReadSet {
	var ret ReadSet
	var numUnique int
	var total int

	for _, r := range r.Reads {
		otherRead := others.GetByName(r.Name)
		if otherRead == nil || r.Score > otherRead.Score {
			ret.Reads = append(ret.Reads, r)

			if verbose {
				if otherRead == nil {
					fmt.Printf("%s is unique. Score=%d\n", r.Name, r.Score)
					numUnique++
				} else {
					difference := r.Score - otherRead.Score
					fmt.Printf("%s is better. %d > %d by %d\n", r.Name,
						r.Score, otherRead.Score, difference)
					total += difference
				}
			}
		}
	}
	ret.BuildNameIndex()

	if verbose {
		fmt.Printf("Unique: %d. Total Score Improvement: %d\n",
			numUnique, total)
	}

	return &ret
}

func ParseReads(fname string, minScore int) *ReadSet {
	var ret ReadSet
	ret.Name = fname
	ret.Reads = make([]Read, 0)

	utils.Lines(fname, func(line string, err error) bool {
		fields := strings.Split(line, "\t")
		if len(fields) < 12 {
			return true
		}
		name := fields[0]
		score := utils.Atoi(fields[11][5:])
		if score >= minScore {
			ret.Reads = append(ret.Reads, Read{name, score})
		}
		return true
	})

	ret.BuildNameIndex()
	return &ret
}

func main() {
	var (
		minScore int
		verbose  bool
	)

	flag.IntVar(&minScore, "min", -40, "Minimum score")
	flag.BoolVar(&verbose, "v", false, "verbose")
	flag.Parse()

	args := flag.Args()

	if len(args) != 2 {
		fmt.Printf("Need two sam files")
		flag.PrintDefaults()
		return
	}

	aN, bN := args[0], args[1]
	a := ParseReads(aN, minScore)
	b := ParseReads(bN, minScore)

	aBetter := a.Better(b, verbose)
	fmt.Printf("%d are better in %s than %s\n", len(aBetter.Reads), aN, bN)

	bBetter := b.Better(a, verbose)
	fmt.Printf("%d are better in %s than %s\n", len(bBetter.Reads), bN, aN)
}
