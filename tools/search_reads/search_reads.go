package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"os"
)

func ParseFastq(fname string, cb func(readData []byte)) {
	var interested bool
	utils.Lines(fname, func(line string, lineError error) bool {
		if line[0] == '@' {
			interested = true
			return true
		}
		if interested {
			cb([]byte(line))
			interested = false
		}
		return true
	})
}

func showMatch(a, b []byte) {
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
		tol     float64
		patS    string
		verbose bool
	)

	flag.BoolVar(&verbose, "v", false, "Verbose")
	flag.Float64Var(&tol, "tol", 0.1, "Tolerance")
	flag.StringVar(&patS, "p", "", "Pattern to look for")
	flag.Parse()

	if len(flag.Args()) != 1 {
		flag.PrintDefaults()
		return
	}

	pattern := []byte(patS)
	n := len(pattern)
	matches := 0

	ParseFastq(flag.Arg(0), func(data []byte) {
		var g genomes.Genomes
		g.Nts = [][]byte{data}
		for search := genomes.NewLinearSearch(&g,
			0, pattern, tol); !search.End(); search.Next() {
			pos, _ := search.Get()
			if verbose {
				showMatch(g.Nts[0][pos:pos+n], pattern)
			}
			matches++
		}
	})
	if verbose {
		fmt.Println(matches)
	}
	if matches == 0 {
		os.Exit(-1)
	}
}
