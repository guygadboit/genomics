package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"os"
)

func ParseFastq(fname string, cb func(name string, readData []byte)) {
	var interested bool
	var name string
	utils.Lines(fname, func(line string, lineError error) bool {
		if line[0] == '@' {
			interested = true
			name = line
			return true
		}
		if interested {
			cb(name, []byte(line))
			interested = false
		}
		return true
	})
}

func showMatch(name string, a, b []byte) {
	fmt.Println(name)
	fmt.Println(string(a))
	fmt.Println(string(b))
	for i := 0; i < len(a); i++ {
		if a[i] == b[i] {
			fmt.Printf("*")
		} else {
			fmt.Printf(" ")
		}
	}
	fmt.Printf("\n\n")
}

func main() {
	var (
		tol        float64
		patS       string
		verbose    bool
		minMatches int
	)

	flag.BoolVar(&verbose, "v", false, "Verbose")
	flag.Float64Var(&tol, "tol", 0.1, "Tolerance")
	flag.StringVar(&patS, "p", "", "Pattern to look for")
	flag.IntVar(&minMatches, "m", 1, "Minimum number of matches")
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.PrintDefaults()
		return
	}

	pattern := []byte(patS)
	n := len(pattern)
	matches := 0

	for _, fname := range flag.Args() {
		ParseFastq(fname, func(name string, data []byte) {
			var g genomes.Genomes
			g.Nts = [][]byte{data}
			for search := genomes.NewLinearSearch(&g,
				0, pattern, tol); !search.End(); search.Next() {
				pos, _ := search.Get()
				if verbose {
					showMatch(name, pattern, g.Nts[0][pos:pos+n])
				}
				matches++
			}
		})
	}
	if verbose {
		fmt.Println(matches)
	}
	if matches < minMatches {
		os.Exit(-1)
	}
}
