package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/utils"
)

type Result struct {
	Aa            byte
	Occurrences   int
	SilentMuts    int
	NonSilentMuts int
}

// We're asking whether some AAs tend to mutate more than others.
type Results map[byte]Result

func main() {
	db := database.NewDatabase()
	results := make(Results)

	g := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	for _, c := range genomes.Translate(g, 0) {
		result := results[c.Aa]
		result.Occurrences++
		for i := 0; i < 3; i++ {
			pos := utils.OneBasedPos(c.Pos + i + 1)
			ids := db.MutationIndex[pos]
			for id, _ := range ids {
				record := db.Get(id)
				for _, mut := range record.NucleotideChanges {
					if mut.Pos == pos {
						if mut.Silence == utils.SILENT {
							result.SilentMuts++
						} else {
							result.NonSilentMuts++
						}
					}
				}
			}
		}
		results[c.Aa] = result
	}

	fmt.Println("AA\tcount\tS\tNS\ttotal\trate\n")
	for k, v := range results {
		total := v.SilentMuts + v.NonSilentMuts
		rate := float64(total)/float64(v.Occurrences)
		fmt.Printf("%c\t%d\t%d\t%d\t%d\t%.2f\n",
			k, v.Occurrences, v.SilentMuts, v.NonSilentMuts, total, rate)
	}
}
