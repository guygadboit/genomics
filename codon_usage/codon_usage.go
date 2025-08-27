package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"math"
	"regexp"
	"strings"
)

type CodonFreq struct {
	Codon        string
	CountPer1000 float64
	Aa           byte
	RSCU         float64 // Relative Synonymous Codon Usage

	// The "relative adaptiveness" is like a normalized version of the RSCU.
	// The most preferred codon will end up as 1.0
	RelAd float64
}

type CodonFreqTable struct {
	Frequencies map[string]CodonFreq
	Optima      map[byte]float64 // The best RSCU for each AA
}

// What's the highest RSCU for a particular AA?
func (ct CodonFreqTable) CalcOptimum(aa byte) float64 {
	syns := genomes.ReverseCodonTable[aa]
	var best float64
	for _, syn := range syns {
		best = max(best, ct.Frequencies[syn].RSCU)
	}
	return best
}

func (ct *CodonFreqTable) UpdateOptima() {
	ct.Optima = make(map[byte]float64)
	for k, _ := range genomes.ReverseCodonTable {
		ct.Optima[k] = ct.CalcOptimum(k)
	}
}

// If any codons are missing, give them RSCUs of 0.5 (which is what Sharp/Li
// says to do)
func (ct *CodonFreqTable) FillInBlanks() {
	for k, v := range genomes.CodonTable {
		if _, there := ct.Frequencies[k]; !there {
			ct.Frequencies[k] = CodonFreq{k, 0, v, 0.5, 0}
		}
	}
}

// Update Optima and Aa, RSCU and RelAd for all the CodonFreqs in
// ct.Frequencies
func (ct *CodonFreqTable) UpdateComputedValues() {
	for k, v := range ct.Frequencies {
		aa := genomes.CodonTable[v.Codon]
		syns := genomes.ReverseCodonTable[aa]

		total := 0.0
		for _, syn := range syns {
			total += ct.Frequencies[syn].CountPer1000
		}

		// What count per 1000 would we expect if all codons synonymous for
		// this AA were used equally?
		expected := total / float64(len(syns))

		// Then the RSCU is just what we do see over that expectation
		rscu := v.CountPer1000 / expected
		ct.Frequencies[k] = CodonFreq{v.Codon, v.CountPer1000, aa, rscu, 0}
	}
	ct.FillInBlanks()

	// Once we know all the RSCUs we can work out the relative adaptations
	ct.UpdateOptima()

	// And now we can work out the relative adaptations
	for k, v := range ct.Frequencies {
		newCf := v
		newCf.RelAd = v.RSCU / ct.Optima[v.Aa]
		ct.Frequencies[k] = newCf
	}
}

// Parse a table pasted out of e.g.
// https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=tisspec
func ParseCUTable(fname string) *CodonFreqTable {
	ret := new(CodonFreqTable)
	ret.Frequencies = make(map[string]CodonFreq)
	utils.Lines(fname, func(line string, err error) bool {
		fields := strings.Split(line, "\t")
		for i := 0; i < len(fields); i += 3 {
			codon := strings.Trim(fields[i], " ")
			if codon == "" {
				continue
			}
			count := utils.Atof(strings.Trim(fields[i+1], " "))
			ret.Frequencies[codon] = CodonFreq{codon, count, 0, 0, 0}
		}
		return true
	})
	ret.UpdateComputedValues()
	return ret
}

// Parse a table pasted out of
// https://www.biologicscorp.com/tools/CAICalculator
func ParseBiologicsTable(fname string) *CodonFreqTable {
	ret := new(CodonFreqTable)
	ret.Frequencies = make(map[string]CodonFreq)
	re := regexp.MustCompile(`\(\s+`)
	utils.Lines(fname, func(line string, err error) bool {
		line = string(re.ReplaceAll([]byte(line), []byte("(")))
		fields := strings.Fields(line)
		if len(fields) != 14 {
			return true
		}
		for i := 0; i < len(fields); i += 7 {
			codon := strings.Trim(fields[i], " ")
			count := utils.Atof(strings.Trim(fields[i+5], " "))
			ret.Frequencies[codon] = CodonFreq{codon, count, 0, 0, 0}
		}
		return true
	})
	ret.UpdateComputedValues()
	return ret
}

func (c *CodonFreqTable) ShowRSCUs() {
	for _, f := range c.Frequencies {
		fmt.Printf("%s %f %f\n", f.Codon, f.CountPer1000, f.RSCU)
	}
}

type RelCodon struct {
	genomes.Codon
	RelAd float64 // The relative adaptation of this codon
}

// A translation decorated with relative adaptation values
type RelTranslation []RelCodon

type RollingAverage struct {
	Total, Count float64
}

// These are defined like this (rather than as pointer methods) because we're
// using them in a map
func (r RollingAverage) Add(val float64) RollingAverage {
	return RollingAverage{r.Total + val, r.Count + 1}
}

func (r RollingAverage) Mean() float64 {
	return r.Total / r.Count
}

func CountSFTypes() {
	counts := make(map[int]int)
	for k, syns := range genomes.ReverseCodonTable {
		if len(syns) == 3 {
			fmt.Printf("%c has 3\n", k)
		}
		counts[len(syns)] += 1
	}
	for k, v := range counts {
		fmt.Printf("%d: %d\n", k, v)
	}
}

// Effective Number of Codons per Wright 1989. This should be a number between
// 20 and 61. 20 means you're very biased about which codons you like, 61 you
// use them evenly.
func (r RelTranslation) ENc() float64 {
	freqs := make(map[string]float64)
	for _, codon := range r {
		freqs[codon.Nts]++
	}

	averageFs := make(map[int]RollingAverage)
	sfTypeCounts := make(map[int]int)

	// Now we consider each AA at a time.
	for aa, syns := range genomes.ReverseCodonTable {
		if aa == '*' { // exclude stop codons
			continue
		}

		// First find n, the total number of times this AA appears
		var n float64
		for _, c := range syns {
			n += float64(freqs[c])
		}
		if n == 0 {
			continue
		}

		// Track how many AAs of each SF type we actually saw to use for the
		// final average calculation below.
		sfType := len(syns)

		// If there's only one synonym we don't count anything because there's
		// nothing to count. But we do still record how many AAs of this SF
		// type we saw
		if sfType == 1 {
			sfTypeCounts[sfType] += 1
			continue
		}

		// Now find the total squared frequency of the synonyms for this AA
		var sum float64
		for _, c := range syns {
			// p is the frequency of this synonym
			p := freqs[c] / n
			sum += p * p
		}

		// Now find the "Homozygosity"
		numerator := (n*sum - 1)
		denom := n - 1
		if numerator == 0 || denom == 0 {
			continue
		}
		sfTypeCounts[sfType]++

		F := numerator / denom

		// And we will need an average for each "SF type" (the "SF type" of a
		// codon is defined by the number of synonyms it has)
		averageFs[sfType] = averageFs[sfType].Add(F)
	}

	// The paper says that if there is no Isoleucine (which is the sole member
	// of sfType 3) then use the average of types 2 and 4 again.
	if _, there := averageFs[3]; !there {
		averageFs[3] = RollingAverage{
			averageFs[2].Total + averageFs[4].Total,
			averageFs[2].Count + averageFs[4].Count,
		}
	}

	// We weight each SF type according to how often it appears in the codon
	// table (not in our sample). This is what the EMBOSS code does (although
	// it isn't clear from the paper)
	weights := map[int]float64{
		2: 9,
		3: 1,
		4: 5,
		6: 3,
	}

	ret := 2.0
	for typ, _ := range sfTypeCounts {
		ra, there := averageFs[typ]
		if !there {
			continue
		}
		ret += float64(weights[typ]) / ra.Mean()
	}

	return min(ret, 61)
}

func CopyWithoutGaps(g *genomes.Genomes, which int) *genomes.Genomes {
	g2 := g.Filter(which)
	g2.DeepCopy(0)
	g2.RemoveGaps()
	return g2
}

/*
Returns what Wright calls "RF1"-- the number of CpGs over the total number
of nucleotides
*/
func CpG(g *genomes.Genomes, which int) float64 {
	// Work with the actual nts ignoring any gaps in the alignment
	g2 := CopyWithoutGaps(g, which)
	nts := g2.Nts[0]

	var count int
	for i := 1; i < len(nts); i++ {
		if nts[i-1] == 'C' && nts[i] == 'G' {
			count++
		}
	}
	return float64(count) / float64(len(nts))
}

// Call cb with a series of the CpG so far at nt position int into the genome
func CumulativeCpG(g *genomes.Genomes, which int, cb func(int, int)) {
	g2 := CopyWithoutGaps(g, which)
	nts := g2.Nts[0]

	var count int
	for i := 1; i < len(nts); i++ {
		if nts[i-1] == 'C' && nts[i] == 'G' {
			count++
		}
		if count > 0 {
			cb(i, count)
		}
	}
}

func calcCAI(trans RelTranslation, ref *CodonFreqTable) float64 {
	prodRSCU, prodBest := 1.0, 1.0
	exceptions := 0

	for _, codon := range trans {
		cf := ref.Frequencies[codon.Nts]
		best := ref.Optima[codon.Aa]

		if cf.RSCU == 0.0 {
			exceptions++
			continue
		}

		prodRSCU += math.Log(cf.RSCU)
		prodBest += math.Log(best)
		// fmt.Println(codon, cf.RSCU, best, prodRSCU, prodBest)

		switch codon.Nts {
		case "ATG":
			fallthrough
		case "TGG":
			exceptions++
		}
	}

	l := float64(len(trans) - exceptions)
	if l == 0 {
		return 0.5
	}
	caiObs := math.Exp(prodRSCU / l)
	caiMax := math.Exp(prodBest / l)

	return caiObs / caiMax
}

func MakeRelTranslation(g *genomes.Genomes,
	which int, ref *CodonFreqTable) (RelTranslation, float64) {
	trans := genomes.Translate(g, which)
	ret := make(RelTranslation, len(trans))

	for i, codon := range trans {
		cf := ref.Frequencies[codon.Nts]
		ret[i] = RelCodon{codon, cf.RelAd}
	}

	return ret, calcCAI(ret, ref)
}

func (r RelTranslation) SubseqCAI(start, end int, ref *CodonFreqTable) float64 {
	return calcCAI(r[start:end], ref)
}

func parseRestrict(g *genomes.Genomes, restrict string) (int, int) {
	for _, orf := range g.Orfs {
		if orf.Name == restrict {
			return orf.Start, orf.End
		}
	}
	ints := utils.ParseInts(restrict, ":")
	return ints[0] - 1, ints[1]
}

func makeGraphData(relTrans RelTranslation, ref *CodonFreqTable,
	which int, window int, labels bool) {
	fname := fmt.Sprintf("%d.txt", which)
	fd, fp := utils.WriteFile(fname)
	defer fd.Close()

	for i := 0; i < len(relTrans)-window; i++ {
		cai := relTrans.SubseqCAI(i, i+window, ref)
		if labels && window == 1 {
			fmt.Fprintf(fp, "%c %f\n", relTrans[i].Aa, cai)
		} else {
			fmt.Fprintf(fp, "%f\n", cai)
		}
	}

	fp.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func writeGnuplotFile(which int, dataName string, name string) {
	fname := fmt.Sprintf("plot-%d.gpi", which)
	fd, fp := utils.WriteFile(fname)
	defer fd.Close()

	fmt.Fprintf(fp, "set title \"%s cumulative CpG count\"\n",
		utils.GnuplotEscape(name))
	fmt.Fprintf(fp, "set term png\n")
	fmt.Fprintf(fp, "set output \"%d.png\"\n", which)
	fmt.Fprintf(fp, "plot \"%s\" with lines notitle\n", dataName)

	fp.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func makeCpGGraphData(g *genomes.Genomes, which int) {
	fname := fmt.Sprintf("%d.txt", which)
	fd, fp := utils.WriteFile(fname)
	defer fd.Close()

	CumulativeCpG(g, which, func(pos int, cpgCount int) {
		fmt.Fprintln(fp, pos+1, cpgCount)
	})

	fp.Flush()
	fmt.Printf("Wrote %s\n", fname)

	writeGnuplotFile(which, fname, g.Names[which])
}

func main() {
	var (
		refName      string
		source, orfs string
		restrict     string
		include      string
		biologics    bool
		graph        bool
		labels       bool
		window       int
		graphCpG     bool
	)

	flag.StringVar(&refName, "ref", "./human", "Reference")
	flag.StringVar(&source, "s",
		"../fasta/SARS2-relatives-short-names.fasta", "Source")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "ORFs")
	flag.StringVar(&include, "i", "0", "Which genomes")
	flag.StringVar(&restrict, "restrict", "", "Restrict to region")
	flag.BoolVar(&biologics, "biologics", false, "Format is biologics")
	flag.BoolVar(&graph, "graph", false, "Write graph data")
	flag.BoolVar(&labels,
		"labels", false, "Label AAs in graph data if window==1")
	flag.IntVar(&window, "window", 1, "Window for graph data")
	flag.BoolVar(&graphCpG, "graph-cpg", false, "Cumulative CpG")
	flag.Parse()

	parse := ParseCUTable
	if biologics {
		parse = ParseBiologicsTable
	}

	ref := parse(refName)

	if orfs == "none" {
		orfs = ""
	}
	g := genomes.LoadGenomes(source, orfs, false)

	if restrict != "" {
		start, end := parseRestrict(g, restrict)
		g.Truncate(start, end)
		// g.Save("check", "check.fasta", 0)
	}

	var indices []int

	if include == "all" {
		indices = make([]int, g.NumGenomes())
		for i := 0; i < len(indices); i++ {
			indices[i] = i
		}
	} else {
		indices = utils.ParseInts(include, ",")
	}

	fmt.Printf("%20s\tCAI\t\tENc\tCpG (%%)\n", "name")
	for _, which := range indices {
		relTrans, cai := MakeRelTranslation(g, which, ref)
		enc := relTrans.ENc()
		cpg := CpG(g, which)
		fmt.Printf("%20s\t%f\t%.3f\t%.3f\n", g.Names[which], cai, enc, cpg*100)
		if graph {
			makeGraphData(relTrans, ref, which, window, labels)
		}
		if graphCpG {
			makeCpGGraphData(g, which)
		}
	}
}
