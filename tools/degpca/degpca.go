package main

import (
	"errors"
	"flag"
	"fmt"
	"genomics/degeneracy"
	"genomics/genomes"
	"genomics/stats"
	"genomics/utils"
	"path"
	"strconv"
	"strings"
)

/*
Contains 8 items, the proportion of codons with 4-fold degeneracy that have A,
C, G, and T in the 3rd position, followed by the same thing for 2-fold
*/
type Row []float64

func countClasses(t degeneracy.Translation) Row {
	ret := make(Row, 8)
	var totalFours, totalTwos float64

	for _, c := range t {
		nt := c.Nts[2]
		fold := c.Fold
		if fold == 3 { // Treat 3-fold as if it were 2
			fold = 2
		}

		var i int
		if fold == 4 {
			i = 0
			totalFours++
		} else {
			i = 4
			totalTwos++
		}

		var j int
		switch nt {
		case 'A':
			j = 0
		case 'C':
			j = 1
		case 'G':
			j = 2
		case 'T':
			j = 3
		}

		ret[i+j]++
	}

	for i := 0; i < 4; i++ {
		ret[i] /= totalFours
	}
	for i := 4; i < 8; i++ {
		ret[i] /= totalTwos
	}

	return ret
}

/*
Load the data from the supplement in the Hassanin paper to check we're on the
same page
*/
func LoadHassanin() [][]float64 {
	ret := make([][]float64, 0)
	currentColumn := -1
	currentRow := 0
	utils.Lines("./supp-table4.txt", func(line string, err error) bool {

		newColumn := true
		switch line {
		case "4X-A":
			currentColumn = 0
		case "4X-C":
			currentColumn = 1
		case "4X-G":
			currentColumn = 2
		case "4X-U":
			currentColumn = 3
		case "2X-A":
			currentColumn = 4
		case "2X-C":
			currentColumn = 5
		case "2X-G":
			currentColumn = 6
		case "2X-U":
			currentColumn = 7
		default:
			newColumn = false
		}

		if currentColumn == -1 {
			return true
		}

		if newColumn {
			fmt.Printf("New column because %s\n", line)
			currentRow = 0
			return true
		}

		if len(ret) < currentRow+1 {
			ret = append(ret, make([]float64, 8))
		}

		val, err := strconv.ParseFloat(line, 64)
		if err == nil {
			fmt.Printf("%d,%d <- %f\n", currentRow, currentColumn, val)
			ret[currentRow][currentColumn] = val
			currentRow++
		}
		return true
	})
	return ret
}

type Source struct {
	fasta   string
	orfs    string
	outName string
	rows    int
}

type Sources []Source

func (s *Sources) Set(value string) error {
	values := strings.Split(value, ",")
	if len(values) != 2 {
		return errors.New("Invalid source specification")
	}

	_, fname := path.Split(values[0])
	outName := utils.BaseName(fname) + ".dat"

	src := Source{values[0], values[1], outName, 0}
	*s = append(*s, src)
	return nil
}

func (s *Sources) String() string {
	return ""
}

func main() {
	var (
		sources Sources
		outName string
	)

	flag.Var(&sources, "s", "List of sources (fasta,orf)")
	flag.StringVar(&outName, "o", "output.dat", "Output file")
	flag.Parse()

	hass := LoadHassanin()
	fmt.Println(hass)
	return

	data := make([][]float64, 0)
	for i, s := range sources {
		g := genomes.LoadGenomes(s.fasta, s.orfs, false)
		for j := 0; j < g.NumGenomes(); j++ {
			t := degeneracy.Translate(g, j)
			classes := countClasses(t)
			data = append(data, classes)
		}
		sources[i].rows = g.NumGenomes()
	}

	result := stats.PCA(2, data)
	fmt.Println("Explained variance ratio:", result.VarianceRatio)

	writeFile := func(fname string, start, end int) {
		fd, fp := utils.WriteFile(fname)
		defer fd.Close()

		for _, row := range result.ReducedData[start:end] {
			fmt.Fprintf(fp, "%f %f\n", row[0], row[1])
		}

		fp.Flush()
		fmt.Printf("Wrote %s\n", fname)
	}

	pos := 0
	for _, s := range sources {
		writeFile(s.outName, pos, pos+s.rows)
		pos += s.rows
	}
}
