package main

import (
	"bufio"
	"encoding/gob"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
	"log"
	"os"
)

const FNAME = "expected-hits.gob"

type ExpectedHits struct {
	Its  int
	Hits []int
}

func (e *ExpectedHits) Init(its int, numGenomes int) {
	e.Its = its
	e.Hits = make([]int, numGenomes)
}

func (e *ExpectedHits) Save() {
	fd, err := os.Create(FNAME)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	enc := gob.NewEncoder(fp)
	err = enc.Encode(e)
	if err != nil {
		log.Fatal(err)
	}
	fp.Flush()
}

func (e *ExpectedHits) Load() error {
	fd, err := os.Open(FNAME)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)
	dec := gob.NewDecoder(fp)
	err = dec.Decode(e)
	if err != nil {
		return err
	}
	return nil
}

/*
Find the expected number of "outgroup" hits for each genome. g should be the
complete set of relatives.
*/
func FindExpected(g *genomes.Genomes, nd *mutations.NucDistro) *ExpectedHits {
	fmt.Println("Doing montecarlo to find expected hits...")
	its := 10000
	var ret ExpectedHits
	ret.Init(its, g.NumGenomes())
	ret.Hits[0] = 0

	for i := 1; i < g.NumGenomes(); i++ {
		g2 := g.Filter(0, i)
		ret.Hits[i] = OutgroupMontecarlo(g2, nd, its)
	}
	ret.Save()
	return &ret
}

func GetExpected(g *genomes.Genomes, nd *mutations.NucDistro) *ExpectedHits {
	var e ExpectedHits
	err := e.Load()
	if err == nil {
		return &e
	}
	return FindExpected(g, nd)
}

/*
For each genome, how many sequences match it above some threshold? Mask are the
ids to mask out: exclude those and any muts matching those from the analysis.
The idea is that we exclude the very close relatives to see what *other*
matches exist not explained by similarity to them
*/
func CountSignificant(
	db *database.Database, ids database.IdSet,
	g *genomes.Genomes,
	expectedHits *ExpectedHits,
	minOR float64,
	maxP float64,
	silent bool,
	mask []int) []int {
	ret := make([]int, g.NumGenomes())

	total := len(ids)
	count := 0
	fmt.Printf("Looking at %d sequences\n", total)

	maskGenomes := g.Filter(mask...)
	maskSet := utils.ToSet(mask)

	for id, _ := range ids {
		r := &db.Records[id]
		for i := 0; i < g.NumGenomes(); i++ {
			if maskSet[i] {
				continue
			}
			g2 := g.Filter(0, i)

			matches := OutgroupMatches(g2,
				r.NucleotideChanges, "", false, silent)

			exclude := OutgroupMatches(maskGenomes,
				r.NucleotideChanges, "", false, silent)

			matches, _ = RemoveIntersection(matches, exclude, true)

			n := len(matches)
			if n == 0 {
				continue
			}
			var ct stats.ContingencyTable
			hits := expectedHits.Hits[i]
			ct.Init(n, len(r.NucleotideChanges)-n, hits, expectedHits.Its-hits)
			OR := ct.CalcOR()
			if OR > minOR {
				_, p := ct.FisherExact()
				if p < maxP {
					fmt.Printf("%s: %s\n", g.Names[i], matches.ToString(false))
					ret[i]++
				}
			}
		}
		count++
		if count%100000 == 0 {
			fmt.Printf("%.0f%% done\n", float64(count*100)/float64(total))
		}
	}

	return ret
}
