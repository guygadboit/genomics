package main

import (
	"genomics/genomes"
	"genomics/mutations"
	"genomics/database"
	"genomics/stats"
	"encoding/gob"
	"log"
	"os"
	"bufio"
	"fmt"
)

const FNAME = "expected-hits.gob"

type ExpectedHits struct {
	Its	int
	Hits	[]int
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

// For each genome, how many sequences match it above some threshold?
func CountSignificant(
	db *database.Database, ids database.IdSet,
	g *genomes.Genomes,
	expectedHits *ExpectedHits,
	minOR float64) []int {
	ret := make([]int, g.NumGenomes())

	total := len(ids)
	count := 0
	fmt.Printf("Looking at %d sequences\n", total)

	for id, _ := range ids {
		r := &db.Records[id]
		for i := 0; i < g.NumGenomes(); i++ {
			g2 := g.Filter(0, i)
			matches := OutgroupMatches(g2, r.NucleotideChanges, "", false)
			n := len(matches)
			if n == 0 {
				continue
			}
			var ct stats.ContingencyTable
			hits := expectedHits.Hits[i]
			ct.Init(n, len(r.NucleotideChanges)-n, hits, expectedHits.Its-hits)
			OR := ct.CalcOR()
			if OR > minOR {
				ret[i]++
			}
		}
		count++
		if count % 100000 == 0 {
			fmt.Printf("%.0f%% done\n", float64(count * 100) / float64(total))
		}
	}

	return ret
}
