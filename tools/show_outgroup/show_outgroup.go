package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/outgroup"
	"genomics/utils"
	"os"
	"path"
	"log"
)

func main() {
	var (
		fasta      string
		window     int
		num        int
		all        bool
		start, end int
	)

	flag.StringVar(&fasta, "fasta", "", "Alignment to use")
	flag.IntVar(&window, "window", 50, "Window Size")
	flag.IntVar(&num, "num", 3, "Number of relatives to look for")
	flag.BoolVar(&all, "all", false, "Ignore num and look at basically all")
	flag.Parse()

	if fasta == "" {
		fasta = path.Join(os.Getenv("GOPATH"),
			"src/genomics/fasta/Hassanin.fasta")
	}

	// start:end (one-based) or a single position
	subseq := flag.Args()[0]

	ss := utils.ParseInts(subseq, ":")
	switch len(ss) {
	case 1:
		start, end = ss[0]-1, ss[0]
	case 2:
		start, end = ss[0]-1, ss[1]
	default:
		log.Fatal("Invalid subsequence")
	}

	g := genomes.LoadGenomes(fasta, "", false)
	start = max(0, start)
	end = min(end, g.Length())

	if all {
		num = g.NumGenomes() - 2
	}

	prox := outgroup.FindNumClosest(g, 0,
		start, end-start, window, num, outgroup.KEEP_BEST)

	for _, p := range prox {
		seq := string(g.Nts[p.Which][start:end])
		fmt.Printf("%d %s: %s\n", p.Which, g.Names[p.Which], seq)
	}
}
