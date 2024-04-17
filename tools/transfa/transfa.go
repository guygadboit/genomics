package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
)

type Mode int

const (
	TRANSLATE Mode = iota
	CODONS
	SIDE_BY_SIDE
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

func showCodons(g *genomes.Genomes, w *bufio.Writer) {
	for i := 0; i < g.NumGenomes(); i++ {
		fmt.Fprintf(w, ">%s\n", g.Names[i])
		t := genomes.Translate(g, i)
		for _, codon := range t {
			fmt.Fprintf(w, "%d %s %c\n",
				codon.Pos, string(codon.Nts), codon.Aa)
		}
	}
}

func main() {
	modes := map[string]Mode{
		"translate": TRANSLATE,
		"codons":    CODONS,
		"ss":        SIDE_BY_SIDE,
	}
	var modeName, orfs, outName string

	flag.StringVar(&modeName, "mode", "translate", "Translation mode")
	flag.StringVar(&orfs, "orfs", "", "ORFs file")
	flag.StringVar(&outName, "o", "", "Output filename")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), orfs, false)
	mode, there := modes[modeName]
	if !there {
		log.Fatal("Invalid mode")
	}

	switch mode {
	case SIDE_BY_SIDE:
		g.SaveWithTranslation(outName, nil)
	case TRANSLATE:
		writeFile(outName, g, translate)
	case CODONS:
		writeFile(outName, g, showCodons)
	}

	if outName != "" {
		fmt.Printf("Wrote %s\n", outName)
	}
}
