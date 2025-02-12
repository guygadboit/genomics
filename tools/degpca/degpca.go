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
func LoadHassanin(p *PCA) {
	data := make(map[string][][]float64)

	utils.Lines("./hassanin-data.txt", func(line string, err error) bool {
		row := make([]float64, 8)
		fields := strings.Fields(line)
		name := fields[0]
		for i, f := range fields[1:] {
			row[i] = utils.Atof(f)
		}
		_, there := data[name]
		if !there {
			data[name] = make([][]float64, 0)
		}
		data[name] = append(data[name], row)
		return true
	})

	for k, v := range data {
		p.AddData(v, k)
	}
}

type Source struct {
	fasta string
	orfs  string
	name  string
	rows  int
}

type Sources []Source

func (s *Sources) Set(value string) error {
	values := strings.Split(value, ",")
	if len(values) != 2 {
		return errors.New("Invalid source specification")
	}

	_, fname := path.Split(values[0])
	name := utils.BaseName(fname)

	src := Source{values[0], values[1], name, 0}
	*s = append(*s, src)
	return nil
}

func (s *Sources) String() string {
	return ""
}

type Label struct {
	name     string
	startRow int
	endRow   int
}

type PCA struct {
	data   [][]float64
	labels []Label
	result stats.PCAResult
}

func NewPCA() *PCA {
	var ret PCA
	ret.data = make([][]float64, 0)
	ret.labels = make([]Label, 0)
	return &ret
}

func (p *PCA) Add(g *genomes.Genomes, name string) {
	start := len(p.data)
	for i := 0; i < g.NumGenomes(); i++ {
		t := degeneracy.Translate(g, i)
		row := countClasses(t)
		p.data = append(p.data, row)
	}
	p.labels = append(p.labels, Label{name, start, len(p.data)})
}

func (p *PCA) AddData(rows [][]float64, name string) {
	start := len(p.data)
	p.data = append(p.data, rows...)
	p.labels = append(p.labels, Label{name, start, len(p.data)})
}

func (p *PCA) Reduce() {
	p.result = stats.PCA(2, p.data)
	fmt.Println("Explained variance ratio:", p.result.VarianceRatio)
}

func (p *PCA) WritePlotData() {
	writeFile := func(fname string, start, end int) {
		fd, fp := utils.WriteFile(fname)
		defer fd.Close()

		for _, row := range p.result.ReducedData[start:end] {
			fmt.Fprintf(fp, "%f %f\n", row[0], row[1])
		}

		fp.Flush()
		fmt.Printf("Wrote %s\n", fname)
	}

	fnames := make([]string, len(p.labels))
	for i, label := range p.labels {
		fname := label.name + ".dat"
		writeFile(fname, label.startRow, label.endRow)
		fnames[i] = fmt.Sprintf("\"%s\"", fname)
	}

	fd, fp := utils.WriteFile("plot.gpi")
	defer fd.Close()

	fmt.Fprintf(fp, "plot %s\n", strings.Join(fnames, ", "))
	fp.Flush()
	fmt.Printf("Wrote plot.gpi\n")
}

func main() {
	var (
		sources Sources
		outName string
		hass    bool
		spikeOnly	bool
	)

	pca := NewPCA()

	flag.Var(&sources, "s", "List of sources (fasta,orf)")
	flag.StringVar(&outName, "o", "output.dat", "Output file")
	flag.BoolVar(&hass, "hass", false, "Include Hassanin data")
	flag.BoolVar(&spikeOnly, "spike", false, "Spike Only")
	flag.Parse()

	if hass {
		if spikeOnly {
			fmt.Printf("Warning: Hassanin data is not spike only\n")
		}
		LoadHassanin(pca)
	}

	for _, s := range sources {
		g := genomes.LoadGenomes(s.fasta, s.orfs, false)
		if spikeOnly {
			S, err := g.Orfs.Find("S")
			if err == nil {
				g.Truncate(S.Start, S.End)
			}
		}
		pca.Add(g, s.name)
	}

	pca.Reduce()
	pca.WritePlotData()
}
