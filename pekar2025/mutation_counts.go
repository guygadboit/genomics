package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/utils"
	"slices"
	"time"
	"bufio"
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

func FindPossibleSilentCT() LocationMap {
	ret := make(LocationMap)
	g := genomes.LoadGenomes("../fasta/WH1.fasta", "../fasta/WH1.orfs", false)
	for _, m := range mutations.PossibleSilentMuts(g, 0) {
		if m.From == 'C' && m.To == 'T' {
			ret[utils.OneBasedPos(m.Pos+1)] = 0
		}
	}
	return ret
}

/*
Make a map of how many of each possible silent C->T mut we see in seqs with
fewer than maxMuts (-1 means unlimited). If after, only look at sequences after
that date. Also return the number of locations with zero sequences having the
C->t mutation there.
*/
func CTDistribution(db *database.Database,
	maxMuts int, after *time.Time) (LocationMap, int) {
	ret := FindPossibleSilentCT()

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
			if m.From == 'C' && m.To == 'T' {
				ret[m.Pos]++
			}

			/*
				if m.Pos == 8782 {
					fmt.Println(r.Summary())
				}
			*/
		}
		return false
	})

	zeros := 0
	for _, v := range ret {
		if v == 0 {
			zeros++
		}
	}
	return ret, zeros
}

func makeGnuplot(maxMuts int, cutoff *time.Time) {
	fd, fp := utils.WriteFile("plot.gpi")
	defer fd.Close()

	var dateS string
	if cutoff != nil {
		dateS = (*cutoff).Format(time.DateOnly)
	} else {
		dateS = "2020-01-01"
	}

	fmt.Fprintf(fp, "set title \"Frequency of each possible C->T silent mut; " +
		"max muts: %d; %s to 2020-12-31\"\n", maxMuts, dateS)

	fmt.Fprintf(fp, "set xlabel \"Nucleotide position\"\n")
	fmt.Fprintf(fp, "set ylabel \"Number of sequences\"\n")

	fmt.Fprintf(fp, "set arrow from 8782, graph -0.1 to 8782, "+
		"graph 0.0 filled lc \"red\"\n")

	fmt.Fprintf(fp, "plot \"output.dat\" with impulses notitle\n")

	fp.Flush()
	fmt.Println("Wrote plot.gpi")
}

func main() {
	var (
		maxMuts  int
		dateS    string
	)

	flag.IntVar(&maxMuts, "m", 1, "maximum number of mutations")
	flag.StringVar(&dateS, "d", "", "only show after date")
	flag.Parse()

	var date *time.Time
	if dateS != "" {
		d := utils.ParseDate(dateS)
		date = &d
	}

	db := database.NewDatabase()
	distro, zeros := CTDistribution(db, maxMuts, date)
	fmt.Printf("%d zeros\n", zeros)

	fd, fp := utils.WriteFile("output.dat")
	defer fd.Close()
	distro.PrintSorted(fp)
	fp.Flush()

	fmt.Println("Wrote output.dat")

	makeGnuplot(maxMuts, date)
}
