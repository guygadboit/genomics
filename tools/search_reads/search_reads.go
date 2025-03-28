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

func printHeading(heading string, char byte) {
	fmt.Println(heading)
	for range heading {
		fmt.Printf("%c", char)
	}
	fmt.Printf("\n")
}

func showMatch(name string, pattern []byte, read []byte, pos int) {

	n := len(pattern)
	a := read[pos : pos+n]

	fmt.Println(name)
	fmt.Println(string(read))
	fmt.Println()

	printHeading("Alignment of match", '=')
	fmt.Println("Read    ", string(a))
	fmt.Println("Pattern ", string(pattern))

	fmt.Printf("         ")
	for i := 0; i < len(a); i++ {
		if a[i] == pattern[i] {
			fmt.Printf("*")
		} else {
			fmt.Printf(" ")
		}
	}
	fmt.Printf("\n\n")
}

/*
pos is the position in read where we found the pattern. refPos is where it is
in the reference. Return the similarity.
*/
func showFullAlignment(read []byte, pos int,
	ref *genomes.Genomes, refPos int) float64 {

	printHeading("Alignment of whole read", '=')
	var lineA, lineB, stars string
	n := len(read)
	matches := 0

	for i, a := range read {
		var b byte
		j := (refPos - pos) + i
		if j >= 0 && j < ref.Length() {
			b = ref.Nts[0][j]
		} else {
			b = '-'
		}
		lineA += string(a)
		lineB += string(b)
		if a == b {
			stars += "*"
			matches++
		} else {
			stars += " "
		}
	}

	for i := 0; i < n; i += 60 {
		end := i + 60
		if end > n {
			end = n
		}
		fmt.Println("Read ", lineA[i:end])
		fmt.Println("Ref  ", lineB[i:end])
		fmt.Println("     ", stars[i:end])
		fmt.Println()
	}
	ss := float64(matches) / float64(n)
	fmt.Printf("Similarity: %.2f\n", ss)
	return ss
}

func main() {
	var (
		tol         float64
		patS        string
		verbose     bool
		minMatches  int
		refName     string
		subseqRange string
	)

	flag.BoolVar(&verbose, "v", false, "Verbose")
	flag.Float64Var(&tol, "tol", 0.1, "Tolerance")
	flag.StringVar(&patS, "p", "", "Pattern to look for")
	flag.StringVar(&refName, "ref", "", "Reference genome")
	flag.StringVar(&subseqRange, "r", "", "Where to find pattern")
	flag.IntVar(&minMatches, "m", 1, "Minimum number of matches")
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.PrintDefaults()
		return
	}

	var pattern []byte

	var start, end int
	if subseqRange != "" {
		ints := utils.ParseInts(subseqRange, ":")
		start = ints[0] - 1
		end = ints[1]
	}

	var ref *genomes.Genomes
	if refName != "" {
		ref = genomes.LoadGenomes(refName, "", false)
		pattern = ref.Nts[0][start:end]
	} else {
		pattern = []byte(patS)
	}

	matches := 0
	for _, fname := range flag.Args() {
		ParseFastq(fname, func(name string, data []byte) {
			var g genomes.Genomes
			g.Nts = [][]byte{data}
			for search := genomes.NewLinearSearch(&g,
				0, pattern, tol); !search.End(); search.Next() {
				pos, _ := search.Get()
				if verbose {
					showMatch(name, pattern, g.Nts[0], pos)
					if ref != nil {
						showFullAlignment(g.Nts[0], pos, ref, start)
					}
				}
				matches++
			}
		})
	}
	if matches < minMatches {
		os.Exit(-1)
	}
	if verbose {
		fmt.Printf("%d reads match\n\n", matches)
	}
}
