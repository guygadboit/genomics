package main

import (
	"fmt"
	"log"
	"strings"
	"genomics/genomes"
	"genomics/align"
	"genomics/database"
)

func convertName(name string) string {
	for _, c := range strings.Split(name, "|") {
		if strings.HasPrefix(c, "EPI_ISL") {
			return c
		}
	}
	return "NOT FOUND"
}

func main() {
	gs := genomes.LoadUnaligned(
		"/fs/f/tmp/GISAID/gisaid_hcov-19_2025_10_31_07.fasta.gz", "", false)
	ref := genomes.LoadGenomes("../../fasta/WH1.fasta",
		"../../fasta/WH1.orfs", false)

	alignment := make([]*genomes.Genomes, 2)
	alignment[0] = ref

	for _, g := range gs {
		alignment[1] = g
		aligned, err := align.Align(alignment, "/fs/f/tmp")
		if err != nil {
			log.Println(err)
			continue
		}

		aligned.Names[1] = convertName(aligned.Names[1])
		r := database.RecordFromAlignment(aligned, 1, nil)

		if len(r.Insertions) != 0 {
			fmt.Println(r.GisaidAccession,
				len(r.Insertions),
				r.Insertions[0].Pos,
				len(r.Insertions[0].Sequence),
				string(r.Insertions[0].Sequence))
		}
	}
}
