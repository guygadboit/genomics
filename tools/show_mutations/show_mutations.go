package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"flag"
	"os"
)

func parseMutation(s string) (pos int, a, b byte) {
	n := len(s)
	pos = utils.Atoi(s[1:n-1]) - 1
	a = s[0]
	b = s[n-1]
	return
}

func main() {
	var fasta string
	var orfs string

	flag.StringVar(&fasta, "f", "", "Fasta file")
	flag.StringVar(&orfs, "orfs", "", "ORFs file")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	for _, mut := range flag.Args() {
		pos, a, b := parseMutation(mut)

		if g.Nts[0][pos] != a {
			fmt.Fprintf(os.Stderr, "%s is invalid\n", mut)
			continue
		}
		err, silent, before, after := genomes.ProteinChange(g, pos,
			0, 0, []byte{b})
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			continue
		}
		if silent {
			fmt.Printf("%s: silent\n", mut)
		} else {
			orfIndex, orfPos, err := g.Orfs.GetOrfRelative(pos)
			if err != nil {
				fmt.Fprintln(os.Stderr, err)
				continue
			}
			name := g.Orfs[orfIndex].Name
			fmt.Printf("%s: %s:%c%d%c\n", mut, name,
				before[0], (orfPos/3)+1, after[0])
		}
	}
}
