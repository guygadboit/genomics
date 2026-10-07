package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/outgroup"
	"genomics/utils"
	"log"
	"os"
	"path"
)

func main() {
	var (
		fasta      string
		window     int
		num        int
		all        bool
		start, end int
		recCA      bool
		xre        bool
	)

	flag.StringVar(&fasta, "fasta", "", "Alignment to use")
	flag.IntVar(&window, "window", 50, "Window Size")
	flag.IntVar(&num, "num", 3, "Number of relatives to look for")
	flag.BoolVar(&all, "all", false, "Ignore num and look at basically all")
	flag.BoolVar(&recCA, "rec", false, "Reconstruct a \"recCA\"")
	flag.BoolVar(&xre, "xre", false, "Exclude BsaI/BsmBI when "+
		"making \"recCA\"")
	flag.Parse()

	if fasta == "" {
		fasta = path.Join(os.Getenv("GOPATH"),
			"src/genomics/fasta/Hassanin.fasta")
	}

	g := genomes.LoadGenomes(fasta, "", false)
	if all {
		num = g.NumGenomes() - 2
	}

	if recCA {
		recCA := outgroup.FindRecCA(g, 0, 1, window, num, xre)
		recCA.SaveMulti("recCA.fasta")
		fmt.Println("Wrote recCA.fasta")
	}

	for _, subseq := range flag.Args() {

		// start:end (one-based) or a single position
		ss := utils.ParseInts(subseq, ":")
		switch len(ss) {
		case 1:
			start, end = ss[0]-1, ss[0]
		case 2:
			start, end = ss[0]-1, ss[1]
		default:
			log.Fatal("Invalid subsequence")
		}

		start = max(0, start)
		end = min(end, g.Length())

		prox := outgroup.FindNumClosest(g, 0,
			start, end-start, window, num, outgroup.KEEP_BEST)

		fmt.Printf("%s:%s\n", subseq, string(g.Nts[0][start:end]))
		for _, p := range prox {
			seq := string(g.Nts[p.Which][start:end])
			fmt.Printf("%d %s: %s (%d/%d)\n",
				p.Which, g.Names[p.Which], seq, p.Differences, p.Window)
		}
	}
}
