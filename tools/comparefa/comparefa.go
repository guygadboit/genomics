package main

import (
	"bufio"
	"flag"
	"fmt"
	. "genomics/comparison"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
	"os/exec"
)

func RunGnuplot(fname string, g *genomes.Genomes, a, b int) {
	gpName := "comparefa-plot.gpi"
	f, err := os.Create(gpName)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	/*
		We don't plot the ratio by default, but the data is there in the file in
		column 5. So if you want to plot that, just edit the gpi script afterwards
	*/

	cleanName := func(s string) string {
		return utils.GnuplotEscape(utils.Shorten(s, 16))
	}

	fmt.Fprintf(w, `
set title "%s vs %s"
set xlabel "nucleotide offset"
set ylabel "count"

plot "cumulative-muts.txt" using 1 title "silent" with lines, \
	"cumulative-muts.txt" using 2 title "non-silent" with lines, \
	"cumulative-muts.txt" using 3 title "insertions" with lines, \
	"cumulative-muts.txt" using 4 title "deletions" with lines
`, cleanName(g.Names[a]), cleanName(g.Names[b]))
	w.Flush()

	cmd := exec.Command("gnuplot", "--persist", gpName)
	err = cmd.Run()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
}

func parseRestrict(g *genomes.Genomes, restrict string) (int, int) {
    for _, orf := range g.Orfs {
        if orf.Name == restrict {
            return orf.Start, orf.End
        }
    }
    ints := utils.ParseInts(restrict, ":")
    return ints[0]-1, ints[1]
}

func main() {
	var (
		orfName         string
		include         string
		saveTrans       string
		oneLine         bool
		keepGaps        bool
		graphData       bool
		silentOnly      bool
		protein         bool
		highlightString string
		highlightFile   string
		showIndels      bool
		showTransitions bool
        restrict        string
	)

	flag.StringVar(&orfName, "orfs", "", "ORFs")
	flag.StringVar(&include, "i", "", "Genomes to include (unset means all)")
	flag.StringVar(&saveTrans, "s", "", "Save translation to file")
	flag.BoolVar(&oneLine, "1", false, "One line output")
	flag.BoolVar(&keepGaps, "gaps", false, "Keep gaps in first genome")
	flag.BoolVar(&graphData, "g", false, "Graph data")
	flag.BoolVar(&silentOnly, "silent", false, "Silent only")
	flag.BoolVar(&protein, "p", false, "input is already translated")
	flag.StringVar(&highlightString, "highlights",
		"", "1-based positions to highlight separated with ,")
	flag.StringVar(&highlightFile, "highlight-file",
		"", "1-based positions to highlight one per line in a file")
	flag.BoolVar(&showIndels, "indels", false, "Show indels")
	flag.BoolVar(&showTransitions, "trans", false, "Show transition counts")
    flag.StringVar(&restrict, "restrict", "", "Restrict to range")
	flag.Parse()

	var g *genomes.Genomes
	for i, fname := range flag.Args() {
		if i == 0 {
			g = genomes.LoadGenomes(fname, orfName, false)
			if !keepGaps {
				g.RemoveGaps()
			}
		} else {
			g2 := genomes.LoadGenomes(fname, orfName, false)
			err := g.AlignCombine(g2)
			if err != nil {
				log.Print(err)
			}
		}
	}

	if orfName == "" {
		g.Orfs = []genomes.Orf{genomes.Orf{0, g.Length(), ""}}
	}

	err := g.CheckOrfs()
	if err != nil {
		fmt.Fprintln(os.Stderr, "Maybe need -gaps?")
		log.Fatal(err)
	}

    if restrict != "" {
        start, end := parseRestrict(g, restrict)
        g.Truncate(start, end)
    }

	var which []int
	if include != "" {
		which = utils.ParseInts(include, ",")
		for _, v := range which {
			if v < 0 || v > g.NumGenomes()-1 {
				log.Fatalf("Invalid index %d", v)
			}
		}
	} else {
		which = make([]int, g.NumGenomes())
		for i := 0; i < g.NumGenomes(); i++ {
			which[i] = i
		}
	}

	for _, w := range which[1:] {
		var c Comparison
		if protein {
			c = CompareProtein(g, which[0], w)
		} else {
			c = Compare(g, which[0], w)
		}
		if oneLine {
			c.OneLineSummary()
		} else if silentOnly {
			c.SilentSummary()
		} else {
			c.Summary(showIndels)
		}
		if len(which) == 2 && graphData {
			fname := "cumulative-muts.txt"
			c.GraphData(fname)
			fmt.Printf("Wrote %s\n", fname)
			RunGnuplot(fname, g, which[0], which[1])
		}
		if showTransitions {
			c.ShowTransitions()
		}
	}

	if len(which) == 2 && saveTrans != "" {
		var highlights []genomes.Highlight
		if highlightFile != "" {
			highlights, err = genomes.ParseHighlightFile(highlightFile,
				true, 'v')
			if err != nil {
				log.Fatal(err)
			}
		} else if highlightString != "" {
			highlights = genomes.ParseHighlights(highlightString,
				",", true, 'v')
		}
		g.SaveWithTranslation(saveTrans, highlights, which[0], which[1])
		fmt.Printf("Wrote %s\n", saveTrans)
	}
}
