package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"strings"
)

type Read struct {
	Name	string
	Pos		int
	Nts		[]byte
}

func (r *Read) Parse(line string) {
	fields := strings.Split(line, "\t")
	r.Name = fields[0]
	r.Pos = utils.Atoi(fields[3]) - 1
	r.Nts = []byte(fields[9])
}

func main() {
	g := genomes.LoadGenomes("../../../fasta/more_relatives.fasta", "", false)

	utils.Lines("472.sam", func(line string, err error) bool {
		var r Read
		if line[0] == '@' {
			return true
		}
		r.Parse(line)
		pos := g.ConvertPosition(472, r.Pos)
		fmt.Printf("pos was %d converted to %d\n", r.Pos, pos)

		var refMatches, otherMatches int
		for i, nt := range r.Nts {
			if g.Nts[0][pos+i] == nt {
				refMatches++
			}
			if g.Nts[472][pos+i] == nt {
				otherMatches++
			}
		}
		fmt.Println(refMatches, otherMatches, len(r.Nts))
		return true
	})
}
