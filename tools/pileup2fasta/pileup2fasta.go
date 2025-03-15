/*
Run this on a pileup generated with e.g.:

$ samtools mpileup sorted-output.sam

And the reference genome that the sam file was originally generated from. It
outputs an alignment of the two in fasta format
*/
package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
	"genomics/utils"
	"log"
	"strings"
)

func showPileup(pu *pileup.Pileup, onlyPos []int) {
	displayRecord := func(pos int) {
		recordI, there := pu.Index[pos]
		if !there {
			fmt.Printf("%d:\n")
			return
		}
		record := &pu.Records[recordI]
		items := make([]string, len(record.Reads))
		for i, read := range record.Reads {
			items[i] = fmt.Sprintf("%cx%d", read.Nt, read.Depth)
		}
		fmt.Printf("%d: %s\n", record.Pos+1, strings.Join(items, ", "))
	}

	if onlyPos != nil {
		for _, pos := range onlyPos {
			displayRecord(pos)
		}
	} else {
		for i := 0; i < pu.MaxPos; i++ {
			displayRecord(i)
		}
	}
}

func main() {
	var (
		reference, output    string
		verbose, veryVerbose bool
		justShow             bool
		showPosS             string
	)

	flag.StringVar(&reference, "ref", "", "Reference genome")
	flag.StringVar(&output, "o", "output.fasta", "Output name")
	flag.BoolVar(&verbose, "v", false, "verbose")
	flag.BoolVar(&veryVerbose, "vv", false, "very verbose")
	flag.BoolVar(&justShow, "show", false, "Just show counts in the pileup")
	flag.StringVar(&showPosS, "pos", "", "Positions to show")
	flag.Parse()

	if len(flag.Args()) != 1 {
		flag.PrintDefaults()
		return
	}

	pileup, err := pileup.Parse(flag.Args()[0])
	if err != nil {
		log.Fatal(err)
	}

	var showPositions []int
	if showPosS != "" {
		justShow = true
		showPositions = utils.ParseInts(showPosS, ",")
		for i, v := range showPositions {
			showPositions[i] = v - 1
		}
	}

	if justShow {
		showPileup(pileup, showPositions)
		return
	}

	g := genomes.LoadGenomes(reference, "", false)

	// The second genome will be what we find in the pileup
	g = g.Filter(0, 0)
	g.DeepCopy(1)
	g.Names[1] = "Pileup"

	nts := g.Nts[1]
	var j int
	for i := 0; i < g.Length(); i++ {
		if j < len(pileup.Records) && pileup.Records[j].Pos == i {
			nts[i] = pileup.Records[j].Reads[0].Nt

			doPrint := veryVerbose || (verbose && nts[i] != g.Nts[0][i])

			if doPrint {
				for _, r := range pileup.Records[j].Reads {
					fmt.Printf("%d%c %d\n", i+1, r.Nt, r.Depth)
				}
			}

			j++
		} else {
			nts[i] = '-'
		}
	}

	if output != "" {
		g.SaveMulti(output)
		fmt.Printf("Wrote %s\n", output)
	}
}
