package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/pileup"
	"genomics/utils"
	"slices"
	"time"
)


type CountCT struct {
	minDepth int
	counts   map[int]int // How many in each position
	silentCT []int       // All possible silent C->T positions
	cutoff   time.Time
	byDate   map[int]int
}

func (c *CountCT) Init(ref *genomes.Genomes, minDepth int, cutoff time.Time) {
	c.minDepth = minDepth
	c.cutoff = cutoff
	muts := mutations.PossibleSilentMuts(ref, 0)
	c.silentCT = make([]int, 0)
	for _, mut := range muts {
		if mut.From == 'C' && mut.To == 'T' {
			c.silentCT = append(c.silentCT, mut.Pos)
		}
	}
	c.counts = make(map[int]int)
	c.byDate = make(map[int]int)
}
func (c *CountCT) Process(record *database.Record, pu *pileup.Pileup) {
	for _, pos := range c.silentCT {
		pur := pu.Get(pos)
		if pur == nil {
			continue
		}
		if !c.cutoff.IsZero() {
			if record.CollectionDate.Compare(c.cutoff) > 0 {
				continue
			}
		}
		for _, read := range pur.Reads {
			if read.Depth >= c.minDepth && read.Nt == 'T' {
				c.counts[pur.Pos] += read.Depth
				delta := int(record.CollectionDate.Sub(utils.Date(2020,
					1, 1)).Hours() / 24)
				c.byDate[delta] = len(c.counts)
			}
		}
	}
}

func (c *CountCT) Display() {
	og := NewOutgroup()

	type count struct {
		pos, count int
	}
	counts := make([]count, 0, len(c.counts))
	for k, v := range c.counts {
		counts = append(counts, count{k, v})
	}
	slices.SortFunc(counts, func(a, b count) int {
		if a.count > b.count {
			return -1
		}
		if a.count < b.count {
			return 1
		}
		return 0
	})
	fmt.Printf("%d possible silent C->Ts\n", len(c.silentCT))
	for _, count := range counts {
		fmt.Printf("%d %d %s\n", count.pos+1,
			count.count, og.Get(count.pos).ToString())
	}

	fmt.Println("By Date")

	type dateCount struct {
		days int // since 2020-01-01
		seen int // number of different C->T muts seen
	}

	dateCounts := make([]dateCount, 0, len(c.byDate))
	for k, v := range c.byDate {
		dateCounts = append(dateCounts, dateCount{k, v})
	}
	slices.SortFunc(dateCounts, func(a, b dateCount) int {
		if a.days < b.days {
			return -1
		}
		if a.days > b.days {
			return 1
		}
		return 0
	})
	for _, dc := range dateCounts {
		fmt.Println(dc.days, dc.seen)
	}
}
