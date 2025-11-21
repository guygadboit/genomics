package main

import (
	"flag"
	"genomics/pileup"
	"log"
)

func main() {
	var reparse bool

	flag.BoolVar(&reparse, "reparse", false, "Parse pileup2fasta "+
		" output rather than an mpileup file")
	flag.Parse()

	var parse func(string) (*pileup.Pileup, error)
	if reparse {
		parse = pileup.Parse2
	} else {
		parse = pileup.Parse
	}

	var puOut *pileup.Pileup

	for i, fname := range flag.Args() {
		pu, err := parse(fname)
		if err != nil {
			log.Fatal(err)
		}
		if i == 0 {
			puOut = pu
		} else {
			puOut = pileup.Combine(puOut, pu)
		}
	}
	puOut.Show(nil)
}
