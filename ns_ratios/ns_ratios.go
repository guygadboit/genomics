package main

import (
	"fmt"
	"genomics/comparison"
	"genomics/genomes"
	"genomics/utils"
	"strings"
)

func LoadShortNames() []string {
	ret := make([]string, 0)
	utils.Lines("../fasta/short_names.txt", func(line string, err error) bool {
		fields := strings.Split(line, " ")
		ret = append(ret, fields[1])
		return true
	})
	return ret
}

var SHORT_NAMES []string

func init() {
	SHORT_NAMES = LoadShortNames()
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	g.Truncate(21562, 25384) // spike only

    /*
    // Interestingly none of these have nearly as high a S/N
	g := genomes.LoadGenomes("../fasta/SARS1-relatives.fasta",
		"../fasta/SARS1.orfs", false)
	g.Truncate(21492, 25259) // spike only (SARS1)
    */

	for i := 0; i < g.NumGenomes(); i++ {
		for j := i + 1; j < g.NumGenomes(); j++ {
			c := comparison.Compare(g, i, j)
			S, NS, _ := c.SilentCount()
			ratio := float64(S) / float64(NS)

			fmt.Printf("%.2f %d %d %s vs %s\n", ratio,
				NS, S, SHORT_NAMES[i], SHORT_NAMES[j])

            /*
			fmt.Printf("%.2f %d %d %s vs %s\n", ratio,
				NS, S, g.Names[i], g.Names[j])
            */
		}
	}
}
