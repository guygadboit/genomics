package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
	"slices"
)

type Mode int

const (
	TRANSLATE Mode = iota
	CODONS
	FREQ
	SIDE_BY_SIDE
	DONT
)

func writeFile(fname string, g *genomes.Genomes,
	cb func(*genomes.Genomes, *bufio.Writer)) {
	var f *os.File
	var err error

	if fname == "" {
		f = os.Stdout
	} else {
		f, err = os.Create(fname)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
	}

	w := bufio.NewWriter(f)
	cb(g, w)
	w.Flush()
}

func translate(g *genomes.Genomes, w *bufio.Writer) {
	for i := 0; i < g.NumGenomes(); i++ {
		if len(g.Nts[i]) < 3 {
			continue
		}
		fmt.Fprintf(w, ">%s\n", g.Names[i])
		t := genomes.Translate(g, i)
		aas := make([]byte, len(t))
		for i, codon := range t {
			aas[i] = codon.Aa
		}
		utils.Wrap(w, aas)
	}
}

type CodonFreq struct {
	nts   string
	count int
}

func showCodons(g *genomes.Genomes, w *bufio.Writer, showFreq bool) {
	counts := make(map[string]int)
	showCodons := !showFreq

	for i := 0; i < g.NumGenomes(); i++ {
		if showCodons {
			fmt.Fprintf(w, ">%s\n", g.Names[i])
		}
		t := genomes.Translate(g, i)
		for _, codon := range t {
			counts[codon.Nts]++
			if showCodons {
				fmt.Fprintf(w, "%d %s %c\n",
					codon.Pos+1, string(codon.Nts), codon.Aa)
			}
		}
	}

	if showFreq {
		fmt.Fprintln(w, "Codon Usage Table")
		for k, v := range genomes.ReverseCodonTable {
			fmt.Fprintf(w, "%c\n", k)
			values := make([]CodonFreq, 0)
			total := 0
			for _, alt := range v {
				c := counts[alt]
				values = append(values, CodonFreq{alt, c})
				total += c
			}
			slices.SortFunc(values, func(a, b CodonFreq) int {
				if a.count < b.count {
					return 1
				}
				if a.count > b.count {
					return -1
				}
				return 0
			})

			for _, v := range values {
				fmt.Fprintf(w, "%s: %d (%.2f)\n",
					v.nts, v.count, float64(v.count)/float64(total))
			}
		}
	}
}

func main() {
	modes := map[string]Mode{
		"translate": TRANSLATE,
		"codons":    CODONS,
		"freq":      FREQ,
		"ss":        SIDE_BY_SIDE,
		"dont":      DONT,
	}
	var (
		modeName, orfs, outName string
		offset                  int
		reverse                 bool
		include                 string
		removeGaps              bool
	)

	flag.StringVar(&modeName, "mode", "translate", "Translation mode")
	flag.StringVar(&orfs, "orfs", "", "ORFs file")
	flag.StringVar(&outName, "o", "", "Output filename")
	flag.IntVar(&offset, "offset", 0, "Frame offset (0, 1 or 2)")
	flag.BoolVar(&reverse, "reverse", false, "Treat as reverse complement")
	flag.StringVar(&include, "i", "", "Genomes to include")
	flag.BoolVar(&removeGaps, "g", true, "Remove gaps from first genome")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), orfs, false)
	if removeGaps {
		g.RemoveGaps()
	}

	if include != "" {
		which := utils.ParseInts(include, ",")
		g = g.Filter(which...)
	}

	mode, there := modes[modeName]
	if !there {
		log.Fatal("Invalid mode")
	}

	if reverse {
		g.Nts[0] = utils.ReverseComplement(g.Nts[0])
	}

	if offset != 0 {
		if offset > 2 {
			log.Fatal("Offset should only be 0, 1 or 2")
		}
		for i, _ := range g.Orfs {
			g.Orfs[i].Start += offset
			g.Orfs[i].End = min(g.Orfs[i].End+offset, g.Length())
		}
	}

	switch mode {
	case SIDE_BY_SIDE:
		fallthrough
	case DONT:
		if outName == "" {
			log.Fatal("Output filename required for ss mode")
		}
	}

	switch mode {
	case SIDE_BY_SIDE:
		g.SaveWithTranslation(outName, nil)
	case TRANSLATE:
		writeFile(outName, g, translate)
	case CODONS:
		fallthrough
	case FREQ:
		writeFile(outName, g, func(g *genomes.Genomes, w *bufio.Writer) {
			showCodons(g, w, mode == FREQ)
		})
	case DONT:
		g.SaveClu(outName, nil)
	}

	if outName != "" {
		fmt.Printf("Wrote %s\n", outName)
	}
}
