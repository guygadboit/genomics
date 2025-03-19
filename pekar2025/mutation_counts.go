package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/utils"
	"log"
	"slices"
	"time"
)

func CTRate(db *database.Database) (int, int) {
	ct, nonCt := 0, 0

	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}

		for _, m := range r.NucleotideChanges {
			if m.From == 'C' && m.To == 'T' {
				ct++
			} else {
				nonCt++
			}
		}
		return false
	})

	fmt.Printf("CT: %d nonCt: %d\n", ct, nonCt)
	return ct, nonCt
}

type LocationMap map[utils.OneBasedPos]int

func (l LocationMap) PrintSorted(fp *bufio.Writer) {
	type result struct {
		pos   utils.OneBasedPos
		count int
	}
	results := make([]result, 0, len(l))

	for k, v := range l {
		results = append(results, result{k, v})
	}

	slices.SortFunc(results, func(a, b result) int {
		if a.count < b.count {
			return 1
		}
		if a.count > b.count {
			return -1
		}
		return 0
	})

	for _, r := range results {
		fmt.Fprintf(fp, "%d %d\n", r.pos, r.count)
	}
}

type Transition struct {
	from, to byte // 'N' means anything
	silent   bool
}

func (t *Transition) Matches(from, to byte) bool {
	fromMatch := t.from == 'N' || t.from == from
	toMatch := t.to == 'N' || t.to == to
	return fromMatch && toMatch
}

func FindPossible(t Transition) LocationMap {
	ret := make(LocationMap)
	g := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)

	if t.silent {
		for _, m := range mutations.PossibleSilentMuts(g, 0) {
			if t.Matches(m.From, m.To) {
				ret[utils.OneBasedPos(m.Pos+1)] = 0
			}
		}
	} else {
		for i := 0; i < g.Length(); i++ {
			if g.Nts[0][i] == t.from {
				ret[utils.OneBasedPos(i+1)] = 0
			}
		}
	}

	return ret
}

/*
Make a map of how many of each possible mut we see in seqs with fewer than
maxMuts (-1 means unlimited). If after, only look at sequences after that date.
Also return the number of locations with zero sequences having the transition
there.
*/
func Distribution(db *database.Database, trans Transition,
	maxMuts int, after *time.Time, verbose bool) (LocationMap, int, int) {
	ret := FindPossible(trans)
	total := len(ret)

	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}

		if maxMuts != -1 && len(r.NucleotideChanges) > maxMuts {
			return false
		}

		if after != nil {
			if r.CollectionDate.Compare(*after) < 0 {
				return false
			}
		}

		for _, m := range r.NucleotideChanges {

			// Disregard if it's not one of the possible silent mutations we
			// already scanned for.
			_, there := ret[m.Pos]
			if !there {
				continue
			}

			if trans.Matches(m.From, m.To) {
				ret[m.Pos]++
				if verbose {
					fmt.Println(r.Summary())
				}
			}
		}
		return false
	})

	zeros := 0
	for _, v := range ret {
		if v == 0 {
			zeros++
		}
	}
	return ret, zeros, total
}

func makeGnuplot(trans Transition, maxMuts int,
	cutoff *time.Time, highlight utils.OneBasedPos) {
	fd, fp := utils.WriteFile("plot.gpi")
	defer fd.Close()

	var dateS string
	if cutoff != nil {
		dateS = (*cutoff).Format(time.DateOnly)
	} else {
		dateS = "2020-01-01"
	}

	var silentS string
	if trans.silent {
		silentS = "silent "
	}

	fmt.Fprintf(fp, "set title \"Frequency of each possible %c->%c "+
		"%smut; max muts: %d; %s to 2020-12-31\"\n",
		trans.from, trans.to, silentS, maxMuts, dateS)

	fmt.Fprintf(fp, "set xlabel \"Nucleotide position\"\n")
	fmt.Fprintf(fp, "set ylabel \"Number of sequences\"\n")

	fmt.Fprintf(fp, "set arrow from %d, graph -0.1 to %d, "+
		"graph 0.0 filled lc \"red\"\n", highlight, highlight)

	fmt.Fprintf(fp, "plot \"output.dat\" with impulses notitle\n")

	fp.Flush()
	fmt.Println("Wrote plot.gpi")
}

func main() {
	var (
		maxMuts    int
		dateS      string
		silent     bool
		transition string
		pos        int
		verbose    bool
	)

	flag.IntVar(&maxMuts, "m", 1, "maximum number of mutations")
	flag.StringVar(&dateS, "d", "", "only show after date")
	flag.BoolVar(&silent, "silent", true, "silent only")
	flag.StringVar(&transition, "t", "CT", "transition to look for")
	flag.IntVar(&pos, "p", 8782, "Position to highlight")
	flag.BoolVar(&verbose, "v", false, "verbose")
	flag.Parse()

	var from, to byte
	if len(transition) != 2 {
		log.Fatal("Invalid transition")
	}
	from, to = transition[0], transition[1]

	highlight := utils.OneBasedPos(pos)

	var date *time.Time
	if dateS != "" {
		d := utils.ParseDate(dateS)
		date = &d
	}

	db := database.NewDatabase()
	trans := Transition{from, to, silent}
	distro, zeros, total := Distribution(db, trans, maxMuts, date, verbose)
	nonZeros := total - zeros
	fmt.Printf("%d/%d zeros (%d non-zeros)\n", zeros, total, nonZeros)

	fd, fp := utils.WriteFile("output.dat")
	defer fd.Close()
	distro.PrintSorted(fp)
	fp.Flush()

	fmt.Println("Wrote output.dat")

	makeGnuplot(trans, maxMuts, date, highlight)
}
