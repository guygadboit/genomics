package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"path"
	"strings"
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
	/*
		fname := "gisaid_hcov-19_2025_10_31_07.fasta.gz"
		fname := "gisaid_hcov-19_2025_10_31_07.fasta.gz"
		fname := gisaid_hcov-19_2025_11_02_09.fasta.gz"
	fname := "gisaid_hcov-19_2025_11_03_08.fasta.gz"
	fname := "2021-08-22.fasta.gz"
	*/
	fname := "B.1.1.7.fasta.gz"
	gs := genomes.LoadUnaligned(path.Join("/fs/f/tmp/GISAID", fname), "", false)
	ref := genomes.LoadGenomes("../../fasta/WH1.fasta",
		"../../fasta/WH1.orfs", false)

	database.RecordsFromUnaligned(gs, ref, nil, convertName,
		func(r *database.Record) {
			for _, ins := range r.Insertions {
				fmt.Println(r.GisaidAccession,
					ins.Pos+1, string(ins.Sequence))
			}
		})

}
