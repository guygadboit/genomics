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
	flag.Parse()

	parse := ParseCUTable
	if biologics {
		parse = ParseBiologicsTable
	}

	ref := parse(refName)
	g := genomes.LoadGenomes(source, orfs, false)

	if restrict != "" {
		start, end := parseRestrict(g, restrict)
		g.Truncate(start, end)
		g.Save("check", "check.fasta", 0)
	}

	for _, which := range utils.ParseInts(include, ",") {
		relTrans, cai := MakeRelTranslation(g, which, ref)
		fmt.Printf("%s: %f\n", g.Names[which], cai)

		if graph {
			makeGraphData(relTrans, ref, which, window, labels)
		}
	}
}
