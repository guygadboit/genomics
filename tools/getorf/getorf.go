package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"log"
)

func makeFname(orfName string, which int, count int) string {
	if count == 1 {
		return fmt.Sprintf("%s.fasta", orfName)
	}
	return fmt.Sprintf("%s-%d.fasta", orfName, which)
}

func main() {
	var (
		removeGaps bool
		orfName    string
	)

	flag.BoolVar(&removeGaps, "g", true, "Remove gaps from first genome")
	flag.StringVar(&orfName, "orf", "S", "The ORF you want")

	flag.Parse()
	args := flag.Args()
	if len(args) < 2 {
		log.Fatal("Need orfs and some fastas")
	}

	for i, arg := range args[1:] {
		g := genomes.LoadGenomes(arg, args[0], false)

		var (
			orf genomes.Orf
			err error
		)

		if orf, err = g.Orfs.Find(orfName); err != nil {
			log.Fatal("Bad ORF name")
		}
		if removeGaps {
			g.RemoveGaps()
		}

		for j := 0; j < g.NumGenomes(); j++ {
			g.Nts[j] = g.Nts[j][orf.Start:orf.End]
		}
		fname := makeFname(orfName, i, len(args[1:]))
		g.SaveMulti(fname)
		fmt.Printf("Wrote %s\n", fname)
	}
}
