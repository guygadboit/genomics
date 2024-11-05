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
	"log"
)

func main() {
	var reference, output string

	flag.StringVar(&reference, "ref", "", "Reference genome")
	flag.StringVar(&output, "o", "output.fasta", "Output name")
	flag.Parse()

	if len(flag.Args()) != 1 {
		flag.PrintDefaults()
		return
	}

	g := genomes.LoadGenomes(reference, "", false)

	pileup, err := pileup.Parse(flag.Args()[0])
	if err != nil {
		log.Fatal(err)
	}

	// The second genome will be what we find in the pileup
	g = g.Filter(0, 0)
	g.DeepCopy(1)
	g.Names[1] = "Pileup"

	nts := g.Nts[1]
	var j int
	for i := 0; i < g.Length(); i++ {
		if j < len(pileup) && pileup[j].Pos == i {
			nts[i] = pileup[j].Nt
			j++
		} else {
			nts[i] = '-'
		}
	}

	g.SaveMulti(output)
	fmt.Printf("Wrote %s\n", output)
}
