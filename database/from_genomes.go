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

func convertAAMuts(c *comparison.Comparison) AAMutations {
	g := c.Genomes
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

func convertInsertions(c *comparison.Comparison) []Insertion {
	g := c.Genomes
	ret := make([]Insertion, 0)
	pos := -1
	seq := make([]byte, 0)
	for _, ins := range c.Insertions {
		nt := g.Nts[c.B][ins]

		// Keep appending to the same insertion if we've got one on the go
		if pos != 1 {
			if ins == pos+len(seq) {
				seq = append(seq, nt)
				continue
			}
		}

		// Otherwise finalize the one we have, if we do have one
		if pos != -1 {
			ret = append(ret, Insertion{pos + 1, seq})
		}

		// And start a new one
		pos = ins
		seq = make([]byte, 1)
		seq[0] = nt
	}

	if pos != -1 {
		ret = append(ret, Insertion{pos + 1, seq})
	}

	return ret
}

func convertDeletions(c *comparison.Comparison) []Range {
	ret := make([]Range, 0)
	current := Range{-1, -1}

	for _, d := range c.Deletions {
		del := utils.OneBasedPos(d + 1)
		// Keep going if we've got a current one on the go
		if del == current.End+1 {
			current.End++
			continue
		}

		// Otherwise finalize the one we have, if we do have one
		if current.Start != -1 {
			ret = append(ret, current)
			current.Start, current.End = -1, -1
		}

		current.Start, current.End = del, del
	}
	if current.Start != -1 {
		ret = append(ret, current)
	}
	return ret
}

func RecordFromAlignment(g *genomes.Genomes, which int,
	getHost func(i int) string) *Record {
	var record Record
	record.GisaidAccession = strings.Split(g.Names[which], " ")[0]

	if getHost != nil {
		record.Host = getHost(which)
	} else {
		record.Host = "ProbablyBat"
	}

	c := comparison.Compare(g, 0, which)

	record.NucleotideChanges = convertNtMuts(&c)
	record.AAChanges = convertAAMuts(&c)
	record.Insertions = convertInsertions(&c)
	record.Deletions = convertDeletions(&c)
	return &record
}

/*
Given an alignment (with WH-1 in the first position) create database records
*/
func AddFromGenomes(db *Database,
	g *genomes.Genomes, getHost func(i int) string) {

	for i := 1; i < g.NumGenomes(); i++ {
		db.Add(RecordFromAlignment(g, i, getHost))
	}
}
