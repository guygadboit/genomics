package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/stats"
	"genomics/utils"
)

/*
Return a contingency table for how well a particular set of mutations matches
g[1:] while not matching mask
*/
func CheckRecord(r *database.Record,
	g *genomes.Genomes, mask *genomes.Genomes, expected MatchOdds,
	silent bool) (stats.ContingencyTable, database.Mutations) {
	var muts database.Mutations
	if silent {
		muts = r.FilterNucleotideChanges(utils.SILENT)
	} else {
		// Silent and non-silent, but we ignore NOT_IN_ORF
		muts = r.FilterNucleotideChanges(utils.SILENT, utils.NON_SILENT)
	}

	matches := OutgroupMatches(g, 1, muts)
	exclude := OutgroupMatches(mask, 0, muts)
	matches, _ = RemoveIntersection(matches, exclude, true)
	n := len(matches)

	var ct stats.ContingencyTable
	var expectedOdds Odds

	if silent {
		expectedOdds = expected.Silent
	} else {
		expectedOdds = expected.All
	}

	ct.Init(n, len(muts)-n, expectedOdds.Matches, expectedOdds.NonMatches)
	ct.CalcOR()
	return ct, matches
}

func CheckCT(ct stats.ContingencyTable,
	r *database.Record, matches database.Mutations,
	name string, minOR, maxP float64) bool {
	OR := ct.OR
	if OR > minOR {
		_, p := ct.FisherExact(stats.GREATER)
		if p < maxP {
			/*
			fmt.Printf("%s %s: %d/%d OR=%.2f p=%.4g\n",
				r.ToString(), name,
				ct.A, ct.A+ct.B, OR, p)
			*/
			fmt.Printf("%s <%s> %s: %s %d/%d OR=%.2f p=%.4g\n",
				r.ToString(), r.SampleSRA, name,
				matches.ToString(),
				ct.A, ct.A+ct.B, OR, p)
			return true
		}
	}
	return false
}

/*
For each genome, how many sequences match it above some threshold? Mask are the
ids to mask out: exclude those and any muts matching those from the analysis.
The idea is that we exclude the very close relatives to see what *other*
matches exist not explained by similarity to them
*/
func CountSignificant(
	db *database.Database, ids []database.Id, g *genomes.Genomes,
	minOR float64, maxP float64, silent bool,
	mask []int, expected []MatchOdds) []int {
	ret := make([]int, g.NumGenomes())

	total := len(ids)
	count := 0
	fmt.Printf("Looking at %d sequences. Silent: %t.\n", total, silent)

	if mask == nil {
		mask = []int{}
	}
	maskGenomes := g.Filter(mask...)
	maskSet := utils.ToSet(mask)

	for _, id := range ids {
		r := &db.Records[id]
		for i := 1; i < g.NumGenomes(); i++ {
			if maskSet[i] {
				continue
			}

			g2 := g.Filter(0, i)
			ct, matches := CheckRecord(r, g2, maskGenomes, expected[i], silent)
			if CheckCT(ct, r, matches, g.Names[i], minOR, maxP) {
				ret[i]++
			}
		}
		count++
		if count%100000 == 0 {
			fmt.Printf("%.0f%% done\n", float64(count*100)/float64(total))
		}
	}

	return ret
}

// Get all the unique mutations in these ids
func GetAllMutations(db *database.Database,
	ids []database.Id) database.Mutations {
	muts := make(map[database.Mutation]bool)
	for _, id := range ids {
		r := &db.Records[id]
		utils.Union(muts, utils.ToSet(r.NucleotideChanges))
	}
	return utils.FromSet(muts)
}

// Like CountSignificant, except look at individual muts rather than whole
// genomes
/*
func CountSignificantMuts(
	db *database.Database, ids []database.Id, g *genomes.Genomes,
	minOR float64, maxP float64, mask []int, expected []MatchOdds) []int {
	ret := make([]int, g.NumGenomes())

	total := len(ids)
	fmt.Printf("Looking at %d sequences\n", total)

	maskGenomes := g.Filter(mask...)
	maskSet := utils.ToSet(mask)
	muts := GetAllMutations(db, ids)

	for i := 1; i < g.NumGenomes(); i++ {
		if maskSet[i] {
			continue
		}
		g2 := g.Filter(0, i)

		// How many of the muts match this genome, without matching any of the
		// ones we're masking out?
		n := 0
		for j, _ := range muts {
			matches := OutgroupMatches(g2, muts[j:j+1], "", false, false)
			exclude := OutgroupMatches(maskGenomes,
				muts[j:j+1], "", false, false)
			matches, _ = RemoveIntersection(matches, exclude, true)
			n += len(matches)
		}
		if n == 0 {
			continue
		}

		var ct stats.ContingencyTable
		expectedOdds := expected[i].All

		ct.Init(n, len(muts)-n, expectedOdds.Matches, expectedOdds.NonMatches)
		OR := ct.CalcOR()
		fmt.Println(i, n, len(muts), OR)
		if OR > minOR {
			_, p := ct.FisherExact()
			if p < maxP {
				ret[i]++
			}
		}
	}

	return ret
}
*/

func FindAllOdds(g *genomes.Genomes, mask []int) []MatchOdds {
	ret := make([]MatchOdds, g.NumGenomes())
	maskGenomes := g.Filter(mask...)
	for i := 1; i < g.NumGenomes(); i++ {
		g2 := g.Filter(0, i)
		ret[i] = OutgroupOdds(g2, maskGenomes)
		fmt.Println(g.Names[i], ret[i])
	}
	return ret
}
