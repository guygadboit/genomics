package main

import (
	"errors"
	"flag"
	"fmt"
	"genomics/degeneracy"
	"genomics/genomes"
	"genomics/stats"
	"genomics/utils"
	"log"
	"path"
	"strings"
)

type Analysis int

const (
	DEG Analysis = iota
	PROT
	NT
	SILENT_NT
	CODON_USAGE_BIAS
)

func (a Analysis) ToString() string {
	switch a {
	case DEG:
		return "Degenerate nt frequencies"
	case PROT:
		return "Protein"
	case NT:
		return "Nucleotide"
	case SILENT_NT:
		return "Silent Nucleotide"
	case CODON_USAGE_BIAS:
		return "Codon usage"
	}
	return "Unknown"
}

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

type SeparateKey string
type SeparateKeys []SeparateKey

func (s *SeparateKeys) Set(value string) error {
	*s = append(*s, SeparateKey(value))
	return nil
}

func (s *SeparateKeys) String() string {
	return ""
}

type Label struct {
	name     string
	startRow int
	endRow   int
}

type PCA struct {
	labels    []Label
	data      [][]float64
	rowLabels []string
	result    stats.PCAResult
}

func NewPCA() *PCA {
	var ret PCA
	ret.data = make([][]float64, 0)
	ret.labels = make([]Label, 0)
	ret.rowLabels = make([]string, 0)
	return &ret
}

func (p *PCA) Add(g *genomes.Genomes, name string) {
	start := len(p.data)
	for i := 0; i < g.NumGenomes(); i++ {
		t := degeneracy.Translate(g, i)
		row := countClasses(t)
		p.data = append(p.data, row)
		p.rowLabels = append(p.rowLabels, fmt.Sprintf("%s (%d)", g.Names[i], i))
	}
	p.labels = append(p.labels, Label{name, start, len(p.data)})
}

func translateAll(g *genomes.Genomes) []genomes.Translation {
	ret := make([]genomes.Translation, g.NumGenomes())
	for i := 0; i < g.NumGenomes(); i++ {
		ret[i] = genomes.Translate(g, i)
	}
	return ret
}

func (p *PCA) makeLabels(g *genomes.Genomes, name string) {
	nGenomes := g.NumGenomes()

	p.rowLabels = make([]string, nGenomes)
	for i := 0; i < nGenomes; i++ {
		p.rowLabels[i] = fmt.Sprintf("%s (%d)", g.Names[i], i)
	}

	labels := make([]Label, 1)
	labels[0] = Label{name, 0, nGenomes}

	p.labels = labels
}

func AAToVector(aa byte) []float64 {
	aas := "ACDEFGHIKLMNPQRSTVWY"
	ret := make([]float64, 20)
	pos := strings.Index(aas, string(aa))
	ret[pos] = 1.0
	return ret
}

func NtToVector(nt byte) []float64 {
	nts := "ACGT"
	ret := make([]float64, 4)
	pos := strings.Index(nts, string(nt))
	ret[pos] = 1.0
	return ret
}

/*
Completely different analysis. Construct a matrix from a single alignment,
where each row contains a vector representing the proteins you have in the
places where any of the genomes differ.
*/
func (p *PCA) AddProtein(g *genomes.Genomes, name string) {
	translations := translateAll(g)
	nGenomes := g.NumGenomes()
	nCodons := len(translations[0])

	data := make([][]float64, nGenomes)

	for i := 0; i < nCodons; i++ {
		here := make(map[byte]bool)
		for j := 0; j < nGenomes; j++ {
			here[translations[j][i].Aa] = true
		}
		if len(here) == 1 {
			continue
		}

		// Ignore insertions/deletions and stop codons
		if here['-'] || here['*'] {
			continue
		}
		for j := 0; j < nGenomes; j++ {
			newCols := AAToVector(translations[j][i].Aa)
			data[j] = append(data[j], newCols...)
		}
	}

	p.data = data
	p.makeLabels(g, name)
}

func (p *PCA) AddNucleotide(g *genomes.Genomes, name string) {
	nGenomes := g.NumGenomes()
	data := make([][]float64, nGenomes)

outer:
	for i := 0; i < g.Length(); i++ {
		here := make(map[byte]bool)
		for j := 0; j < nGenomes; j++ {
			here[g.Nts[j][i]] = true
		}
		if len(here) == 1 {
			continue
		}
		for k, _ := range here {
			if !utils.IsRegularNt(k) {
				continue outer
			}
		}
		for j := 0; j < nGenomes; j++ {
			newCols := NtToVector(g.Nts[j][i])
			data[j] = append(data[j], newCols...)
		}
	}
	p.data = data
	p.makeLabels(g, name)
}

func (p *PCA) AddSilentNucleotide(g *genomes.Genomes, name string) {
	translations := translateAll(g)
	nGenomes := g.NumGenomes()
	nCodons := len(translations[0])

	data := make([][]float64, nGenomes)

	for i := 0; i < nCodons; i++ {
		aasHere := make(map[byte]bool)
		ntsHere := make(map[string]bool)
		for j := 0; j < nGenomes; j++ {
			aasHere[translations[j][i].Aa] = true
			ntsHere[translations[j][i].Nts] = true
		}

		// If the nts are all the same, not interested
		if len(ntsHere) == 1 {
			continue
		}

		// Ignore insertions/deletions and stop codons
		if aasHere['-'] || aasHere['*'] {
			continue
		}

		// But we require the AAs to all be the same, since we're looking at
		// silent nucleotide changes here
		if len(aasHere) != 1 {
			continue
		}

		for j := 0; j < 3; j++ {
			ntsHere := make(map[byte]bool)
			for k := 0; k < nGenomes; k++ {
				ntsHere[translations[k][i].Nts[j]] = true
			}
			if len(ntsHere) == 1 {
				continue
			}
			for k := 0; k < nGenomes; k++ {
				nt := translations[k][i].Nts[j]
				newCols := NtToVector(nt)
				data[k] = append(data[k], newCols...)
			}
		}
	}

	p.data = data
	p.makeLabels(g, name)
}

type CodonFreq struct {
	nts  string
	freq float64
}

/*
Tells you how many times each codon occurs as a ratio of the total number
*/
func codonFreq(t genomes.Translation) []CodonFreq {
	freqs := make(map[string]float64)

	// Make sure we have an entry for each codon, even if it's zero
	for k, _ := range genomes.CodonTable {
		freqs[k] = 0
	}

	for _, codon := range t {
		if _, there := freqs[codon.Nts]; !there {
			// There might be '-' or things in there
			continue
		}
		freqs[codon.Nts]++
	}

	ret := make([]CodonFreq, 0, len(freqs))
	total := float64(len(t))
	for k, v := range freqs {
		ret = append(ret, CodonFreq{k, v / total})
	}

	// Sort them into a canonical order, alphabetical will do
	utils.SortByKey(ret, false, func(cf CodonFreq) string {
		return cf.nts
	})

	return ret
}

func (p *PCA) AddCodons(g *genomes.Genomes, name string) {
	translations := translateAll(g)
	data := make([][]float64, g.NumGenomes())

	for i, t := range translations {
		cfs := codonFreq(t)
		// Now we just use those frequencies, in a canonical order, as our
		// vector.
		v := make([]float64, len(cfs))
		for j, cf := range cfs {
			v[j] = cf.freq
		}
		data[i] = v
	}
	p.data = data
	p.makeLabels(g, name)
}

// Separate out the rows based on a list of names
func (p *PCA) SeparateNames(g *genomes.Genomes, names SeparateKeys) {

	// For each name, which rows belong to that name?
	rows := make(map[string][]int)

	// Rows that didn't match any of the names
	leftover := make([]int, 0)

outer:
	for i := 0; i < g.NumGenomes(); i++ {
		for _, k := range names {
			name := string(k)
			if strings.Contains(g.Names[i], name) {
				if _, there := rows[name]; !there {
					rows[name] = make([]int, 0)
				}
				rows[name] = append(rows[name], i)
				continue outer
			}
		}
		leftover = append(leftover, i)
	}

	// Now rearrange the rows and labels.
	rowLabels := make([]string, 0, g.NumGenomes())
	labels := make([]Label, 0, len(names)+1)
	data := make([][]float64, 0, len(p.data))

	start := 0
	for k, v := range rows {
		for _, index := range v {
			rowLabels = append(rowLabels, p.rowLabels[index])
			data = append(data, p.data[index])
		}
		n := len(v)
		labels = append(labels, Label{k, start, start + n})
		start += n
	}

	for _, index := range leftover {
		rowLabels = append(rowLabels, p.rowLabels[index])
		data = append(data, p.data[index])
	}

	labels = append(labels, Label{"Other", start, len(p.data)})

	p.labels = labels
	p.data = data
	p.rowLabels = rowLabels
}

// Separate out the rows based on a list of indices
func (p *PCA) Separate(g *genomes.Genomes, indices map[int]bool) {
	a := make([][]float64, 0)
	b := make([][]float64, 0)

	aLabels := make([]string, 0)
	bLabels := make([]string, 0)

	for i := 0; i < g.NumGenomes(); i++ {
		if indices[i] {
			a = append(a, p.data[i])
			aLabels = append(aLabels, p.rowLabels[i])
		} else {
			b = append(b, p.data[i])
			bLabels = append(bLabels, p.rowLabels[i])
		}
	}
	name := "special"
	labels := make([]Label, 2)
	labels[0] = Label{name, 0, len(a)}
	labels[1] = Label{fmt.Sprintf("non-%s", name), len(a), g.NumGenomes()}

	data := append(a, b...)
	p.labels = labels
	p.data = data
	p.rowLabels = append(aLabels, bLabels...)
}

func (p *PCA) AddData(rows [][]float64, name string) {
	start := len(p.data)
	p.data = append(p.data, rows...)
	p.labels = append(p.labels, Label{name, start, len(p.data)})

	rowLabels := make([]string, len(rows))
	for i := 0; i < len(rowLabels); i++ {
		rowLabels[i] = fmt.Sprintf("%s-%d\n", name, i)
	}
	p.rowLabels = append(p.rowLabels, rowLabels...)
}

func (p *PCA) Reduce() {
	fmt.Printf("Reducing %dx%d matrix...\n", len(p.data), len(p.data[0]))
	p.result = stats.PCA(2, p.data)
	fmt.Println("Explained variance ratio:", p.result.VarianceRatio)
}

func (p *PCA) WritePlotData(analysis Analysis) {
	writeFile := func(fname string, start, end int) {
		fd, fp := utils.WriteFile(fname)
		defer fd.Close()

		for i, row := range p.result.ReducedData[start:end] {
			fmt.Fprintf(fp, "%f %f # %s\n",
				row[0], row[1], p.rowLabels[start+i])
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

	fmt.Fprintf(fp, "set title \"%s\"\n\n", analysis.ToString())
	fmt.Fprintf(fp, "plot %s\n", strings.Join(fnames, ", "))
	fp.Flush()
	fmt.Printf("Wrote plot.gpi\n")
}

func main() {
	var (
		sources   Sources
		outName   string
		hass      bool
		spikeOnly bool

		sepKeys     SeparateKeys
		exclude     string
		analysisS   string
		sepIndicesS string
		separate    string
		analysis    Analysis
	)

	pca := NewPCA()

	flag.Var(&sources, "s", "List of sources (fasta,orf)")
	flag.StringVar(&outName, "o", "output.dat", "Output file")
	flag.BoolVar(&hass, "hass", false, "Include Hassanin data")
	flag.BoolVar(&spikeOnly, "spike", false, "Spike Only")
	flag.StringVar(&analysisS, "mode", "nt", "deg|prot|nt|snt")
	flag.Var(&sepKeys, "separate", "Strings to separate on (e.g. 'Pangolin')")
	flag.StringVar(&sepIndicesS, "sepint", "", "Indices to separate")
	flag.StringVar(&exclude, "e", "", "Genomes to exclude")
	flag.Parse()

	switch analysisS {
	case "deg":
		analysis = DEG
	case "prot":
		analysis = PROT
	case "nt":
		analysis = NT
	case "snt":
		analysis = SILENT_NT
	case "cub":
		analysis = CODON_USAGE_BIAS
	default:
		log.Fatal("Unrecognized analysis mode")
	}

	prot_or_nt := analysis == PROT || analysis == NT || analysis == SILENT_NT

	if prot_or_nt {
		if len(sources) != 1 {
			log.Fatal("Only use one alignment for this")
		}
	}

	var exIndices map[int]bool

	if exclude != "" {
		if !prot_or_nt {
			log.Fatal("Exclude is currently only for protein/nt PCA")
		}
		exIndices = utils.ToSet(utils.ParseInts(exclude, ","))
	}

	if separate != "" {
		if !prot_or_nt {
			log.Fatal("Separate is only for protein/nt PCA")
		}
	}

	if hass {
		if spikeOnly {
			fmt.Printf("Warning: Hassanin data is not spike only\n")
		}
		LoadHassanin(pca)
	}

	for _, s := range sources {
		g := genomes.LoadGenomes(s.fasta, s.orfs, false)
		g.RemoveGaps()
		if spikeOnly {
			S, err := g.Orfs.Find("S")
			if err == nil {
				g.Truncate(S.Start, S.End)
			}
		}
		if exIndices != nil {
			which := make([]int, 0)
			for i := 0; i < g.NumGenomes(); i++ {
				if !exIndices[i] {
					which = append(which, i)
				}
			}
			g = g.Filter(which...)
		}
		switch analysis {
		case DEG:
			pca.Add(g, s.name)
		case PROT:
			pca.AddProtein(g, s.name)
		case NT:
			pca.AddNucleotide(g, s.name)
		case SILENT_NT:
			pca.AddSilentNucleotide(g, s.name)
		case CODON_USAGE_BIAS:
			pca.AddCodons(g, s.name)
		}
		if sepIndicesS != "" {
			sepIndices := utils.ToSet(utils.ParseInts(sepIndicesS, ","))
			pca.Separate(g, sepIndices)
		} else if len(sepKeys) != 0 {
			pca.SeparateNames(g, sepKeys)
		}
		if separate != "" || sepIndicesS != "" {
			sepIndices := utils.ToSet(utils.ParseInts(sepIndicesS, ","))
			pca.Separate(g, sepIndices)
		}
	}

	pca.Reduce()
	pca.WritePlotData(analysis)
}
