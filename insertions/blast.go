package main

import (
	"fmt"
	"genomics/stats"
	"genomics/utils"
	"log"
	"slices"
	"os"
	"bufio"
	"encoding/gob"
)

type BlastResult struct {
	stats.BlastResult
	id int
}

type BlastResults []BlastResult

// Maps insertion id to a BlastResult. We will do one of these for each
// organism we're interested in.
func BlastInsertions(insertions []Insertion, genome string) []BlastResult {
	bc := stats.BlastDefaultConfig()
	ret := make([]BlastResult, 0)

	for _, ins := range insertions {
		results, _ := stats.Blast(bc, genome,
			ins.Nts, 1, 1, "", stats.NOT_VERBOSE)
		switch len(results) {
		case 0:
			continue
		case 1:
			fmt.Println(ins.Id, results[0].E)
			ret = append(ret, BlastResult{results[0], ins.Id})
		default:
			log.Fatal("Didn't think this should happen")
		}
	}

	slices.SortFunc(ret, func(a, b BlastResult) int {
		if a.E < b.E {
			return -1
		}
		if a.E > b.E {
			return 1
		}
		return 0
	})

	return ret
}


func (br *BlastResults) Save(fname string) {
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	enc := gob.NewEncoder(fp)
	err = enc.Encode(br)
	if err != nil {
		log.Fatal(err)
	}
	fp.Flush()
	fmt.Printf("Blast results saved to %s\n", fname)
}

func (br *BlastResults) Load(fname string) error {
	fd, err := os.Open(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)
	dec := gob.NewDecoder(fp)
	err = dec.Decode(&br)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Blast results loaded from %s\n", fname)
	return nil
}

func ProkBlast(insertions []Insertion,
	filters []filterFunc) BlastResults {
	fname := "blast-results.gob"
	ret := make(BlastResults, 0)
	err := ret.Load(fname)
	if err == nil {
		return ret
	}

	bc := stats.BlastDefaultConfig()
	filterInsertions(insertions, filters, func(ins *Insertion) {
		fmt.Printf("%s: %dnts\n", ins.ToString(), len(ins.Nts))
		results, err := stats.Blast(bc,
			"/fs/f/genomes/blast/prok/ref_prok_rep_genomes",
			ins.Nts, 10, 10, "", stats.NOT_VERBOSE)
		if err != nil {
			log.Print(err)
			return
		}
		for _, r := range results {
			fmt.Println(ins.Id, r.E, r.Organism)
			ret = append(ret, BlastResult{r, ins.Id})
		}

	}, false)

	ret.Save(fname)
	return ret
}

func showResults(results BlastResults) {
	organisms := make(map[string]int)
	for _, r := range results {
		organisms[r.Organism]++
		fmt.Printf("%d %s %d %f\n", r.id, r.Organism, r.Length, r.E)
	}

	type count struct {
		key   string
		value int
	}
	counts := make([]count, 0, len(results))

	for k, v := range organisms {
		counts = append(counts, count{k, v})
	}
	utils.SortByKey(counts, true, func(c count) int {
		return c.value
	})
	for _, c := range counts {
		fmt.Printf("%s: %d\n", c.key, c.value)
	}
}
