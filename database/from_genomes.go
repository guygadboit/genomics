package database

import (
	"genomics/comparison"
	"genomics/genomes"
	"genomics/utils"
	"strings"
)

func convertNtMuts(c *comparison.Comparison) Mutations {
	ret := make(Mutations, len(c.NtMuts))
	for i, mut := range c.NtMuts {
		ret[i] = Mutation{utils.OneBasedPos(mut.Pos + 1),
			mut.A, mut.B, mut.Silence}
	}
	return ret
}

func convertAAMuts(c *comparison.Comparison, g *genomes.Genomes) AAMutations {
	ret := make(AAMutations, len(c.Muts))
	for i, mut := range c.Muts {
		var gene string
		orfI, oPos, err := g.Orfs.GetOrfRelative(mut.Pos)
		if err != nil {
			continue
		}

		gene = g.Orfs[orfI].Name
		ret[i] = AAMutation{
			Mutation{utils.OneBasedPos(oPos + 1),
				mut.A, mut.B, utils.NON_SILENT},
			gene}
	}
	return ret
}

/*
Given an alignment (with WH-1 in the first position) create database records
*/
func AddFromGenomes(db *Database,
	g *genomes.Genomes, getHost func(i int) string) {

	for i := 1; i < g.NumGenomes(); i++ {
		var record Record
		record.GisaidAccession = strings.Split(g.Names[i], " ")[0]

		if getHost != nil {
			record.Host = getHost(i)
		} else {
			record.Host = "ProbablyBat"
		}

		c := comparison.Compare(g, 0, i)

		record.NucleotideChanges = convertNtMuts(&c)
		record.AAChanges = convertAAMuts(&c, g)

		// FIXME: Not doing insertions or deletions for now
		db.Add(&record)
	}
}
