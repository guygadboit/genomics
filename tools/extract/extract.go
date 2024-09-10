package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
	"path/filepath"
)

func makeOutname(fname string, which int, pos int, reverse bool) string {
	_, base := filepath.Split(fname)
	root := utils.BaseName(base)
	var reversed string
	if reverse {
		reversed = "R"
	}
	return fmt.Sprintf("%s-%d-%d%s.fasta", root, which, pos+1, reversed)
}

func writeOutput(fname string, name string, nts []byte) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create output file")
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	fmt.Fprintf(w, ">%s\n", name)
	utils.Wrap(w, nts)
	w.Flush()
}

type LocationSearch struct {
	location  int
	reverse   bool
	end       bool
	genomeLen int
}

func (l *LocationSearch) Start() {}

func (l *LocationSearch) Get() (int, error) {
	return l.location, nil
}

func (l *LocationSearch) IsForwards() bool {
	return !l.reverse
}

func (l *LocationSearch) Next() {
	l.end = true
}

func (l *LocationSearch) End() bool {
	return l.end
}

func (l *LocationSearch) GenomeLength() int {
	return l.genomeLen
}

func main() {
	var (
		bidi      bool
		patString string
		context   int
		location  int
		reverse   bool
		tolerance float64
	)

	flag.StringVar(&patString, "p", "", "Pattern to look for")
	flag.IntVar(&context, "c", 1000, "Number of nt of context")
	flag.BoolVar(&bidi, "b", true, "Look in both directions")
	flag.IntVar(&location, "l", 0, "One-based location (rather than pattern)")
	flag.BoolVar(&reverse, "r", false, "Reverse (if you used location)")
	flag.Float64Var(&tolerance, "tol", 0.0, "Search tolerance")
	flag.Parse()

	pattern := []byte(patString)
	location -= 1

	for _, fname := range flag.Args() {
		g := genomes.LoadGenomes(fname, "", false)
		for i := 0; i < g.NumGenomes(); i++ {

			var search genomes.Search

			if location != -1 {
				search = &LocationSearch{location, reverse, false, g.Length()}
			} else if bidi {
				search = genomes.NewBidiLinearSearch(g, i, pattern, tolerance)
			} else {
				search = genomes.NewLinearSearch(g, i, pattern, tolerance)
			}

			for ; !search.End(); search.Next() {
				pos, _ := search.Get()
				start := pos - context
				if start < 0 {
					start = 0
				}
				end := pos + context + len(pattern)
				if end >= g.Length() {
					end = g.Length()
				}

				reverse := !search.IsForwards()
				subseq := g.Nts[i][start:end]
				if reverse {
					subseq = utils.ReverseComplement(subseq)
				}

				outFname := makeOutname(fname, i, pos, reverse)
				writeOutput(outFname, g.Names[i], subseq)
				fmt.Printf("Wrote %s\n", outFname)
			}
		}
	}
}
