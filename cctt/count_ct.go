package main

import (
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/pileup"
	"genomics/utils"
	"genomics/stats"
	"slices"
	"time"
)

type Contents struct {
	C8782 int // The depth of C at 8782
	T8782 int

	C28144 int
	T28144 int

	Classification string
	Significance   float64
}

func (c *Contents) ToString() string {
	return fmt.Sprintf("%s %d|%d %d|%d %.4g", c.Classification,
		c.C8782, c.T8782, c.C28144, c.T28144, c.Significance)
}

/*
Suppose you find much more C than T in both locations. What is the probability
of doing that under the assumption that what you have is a mix of CT and TC? If
it were the latter, you would expect C/T to roughly equal T/C in the two
locations.
*/
func (c *Contents) CalcSignificance() {
	var ct stats.ContingencyTable
	ct.Init(c.C8782, c.T8782, c.T28144, c.C28144)
	ct.FisherExact(stats.TWO_SIDED)
	c.Significance = ct.P
}

func Classify(pu *pileup.Pileup, minDepth int) Contents {
	pos8782 := pu.Get(8782 - 1)
	pos28144 := pu.Get(28144 - 1)

	ret := Contents{
		pos8782.GetDepthOf('C'), pos8782.GetDepthOf('T'),
		pos28144.GetDepthOf('C'), pos28144.GetDepthOf('T'),
		"-", 1.0,
	}

	if ret.C8782 >= minDepth &&
		ret.T8782 < minDepth &&
		ret.C28144 >= minDepth &&
		ret.T28144 < minDepth {
		ret.Classification = "CC*"
	} else if ret.T8782 >= minDepth &&
		ret.C8782 < minDepth &&
		ret.T28144 >= minDepth &&
		ret.C28144 < minDepth {
		ret.Classification = "TT*"
	} else if ret.C8782 > ret.T8782 &&
		ret.C28144 > ret.T28144 &&
		ret.C8782 >= minDepth &&
		ret.C28144 >= minDepth {
		ret.Classification = "CC>"
	} else if ret.T8782 > ret.C8782 &&
		ret.T28144 > ret.C28144 &&
		ret.T8782 >= minDepth &&
		ret.T28144 >= minDepth {
		ret.Classification = "TT>"
	}
	ret.CalcSignificance()
	return ret
}

func (c *Contents) QSDepth() int {
	return utils.Max(utils.Min(c.C8782, c.T8782), utils.Min(c.C28144, c.T28144))
}



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

func CountReadsCT(db *database.Database,
	ids []database.Id, minDepth int, prefix string, cutoff time.Time) {
	var c CountCT
	ref := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	c.Init(ref, minDepth, cutoff)
	ProcessReads(db, ids, prefix, &c)
	c.Display()
}

