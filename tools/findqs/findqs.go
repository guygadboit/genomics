package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/pileup"
	"genomics/stats"
	"log"
	"os"
	"path"
	"slices"
)

func printSorted(counts map[int]int) {
	type record struct {
		k, v int
	}
	records := make([]record, 0, len(counts))
	total := 0
	for k, v := range counts {
		records = append(records, record{k, v})
		total += v
	}
	slices.SortFunc(records, func(a, b record) int {
		if a.v < b.v {
			return 1
		}
		if a.v > b.v {
			return -1
		}
		return 0
	})

	fmt.Printf("Total matches: %d\n", total)
	for _, rec := range records {
		fmt.Printf("%d: %d\n", rec.k, rec.v)
	}
}

func majority(alleles map[byte]int) byte {
	var ret byte
	best := -1
	for k, v := range alleles {
		if v > best {
			best = v
			ret = k
		}
	}
	return ret
}

func iterateSilent(g *genomes.Genomes, cb func(int, byte, byte)) {
	for _, mut := range mutations.PossibleSilentMuts(g, 0) {
		cb(mut.Pos, mut.From, mut.To)
	}
}

func iterateAll(g *genomes.Genomes, cb func(int, byte, byte)) {
	for i := 0; i < g.Length(); i++ {
		for _, nt := range []byte{'G', 'A', 'T', 'C'} {
			if nt != g.Nts[0][i] {
				cb(i, g.Nts[0][i], nt)
			}
		}
	}
}

func ExpectedMajorityRate(g *genomes.Genomes,
	requireSilent bool, requireTC bool) (int, int, float64) {
	var iterate func(g *genomes.Genomes, cb func(int, byte, byte))
	if requireSilent {
		iterate = iterateSilent
	} else {
		iterate = iterateAll
	}

	var count, total int
	iterate(g, func(pos int, from, to byte) {
		if requireTC {
			if from != 'T' || to != 'C' {
				return
			}
		}
		alleles := make(map[byte]int)
		for i := 1; i < g.NumGenomes(); i++ {
			alleles[g.Nts[i][pos]]++
		}
		if majority(alleles) == to {
			count++
		}
		total++
	})

	fmt.Println(count, total)
	return count, total - count, float64(count) / float64(total)
}

func ExpectedMatchRate(g *genomes.Genomes,
	requireSilent bool, requireTC bool) (int, int, float64) {
	var iterate func(g *genomes.Genomes, cb func(int, byte, byte))
	if requireSilent {
		iterate = iterateSilent
	} else {
		iterate = iterateAll
	}

	var count, total int
	iterate(g, func(pos int, from, to byte) {
		for i := 1; i < g.NumGenomes(); i++ {
			if requireTC {
				if from != 'T' || to != 'C' {
					return
				}
			}
			if g.Nts[i][pos] == to {
				count++
				break
			}
		}
		total++
	})

	fmt.Println(count, total)
	return count, total - count, float64(count) / float64(total)
}

type Outgroup struct {
	totalSilentMuts int
	totalOGMatches  int
	recCA           *genomes.Genomes
}

func (o *Outgroup) Init(recCA, g *genomes.Genomes) {
	o.recCA = recCA

	possible := mutations.PossibleSilentMuts(g, 0)
	o.totalSilentMuts = len(possible)
	for _, mut := range possible {
		if recCA.Nts[0][mut.Pos] == mut.To {
			o.totalOGMatches++
		}
	}
}

func (o *Outgroup) Matches(pos int, nt byte) bool {
	return o.recCA.Nts[0][pos] == nt
}

func (o *Outgroup) IsRemarkable(numMatches, numMuts int) (float64, float64) {
	var ct stats.ContingencyTable
	ct.Init(numMatches, numMuts-numMatches,
		o.totalOGMatches, o.totalSilentMuts-o.totalOGMatches)
	return ct.FisherExact(stats.GREATER)
}

func Compare(pu *pileup.Pileup,
	g *genomes.Genomes, minDepth int,
	minRatio float64, requireSilent bool,
	requireTC bool, showReads bool, verbose bool,
	og *Outgroup) {

	total, totalOGMatches := 0, 0
	rank0Depth := 0

	for i := 0; i < g.Length(); i++ {
		rec := pu.Get(i)
		if rec == nil {
			continue
		}
		for rank, read := range rec.Reads {
			if rank == 0 {
				rank0Depth = read.Depth
			}
			if read.Nt == g.Nts[0][i] {
				continue
			}
			if read.Depth < minDepth {
				break
			}
			if rank > 0 && minRatio != 0.0 {
				if float64(read.Depth)/float64(rank0Depth) < minRatio {
					break
				}
			}
			silent, _, _ := genomes.IsSilentWithReplacement(g,
				i, 0, 0, []byte{read.Nt})
			if requireSilent && !silent {
				continue
			}

			if requireTC {
				if read.Nt != 'C' || g.Nts[0][i] != 'T' {
					continue
				}
			}

			var silentS string
			if silent {
				silentS = "*"
			}

			total++
			matchesOg := og.Matches(i, read.Nt)
			if matchesOg {
				totalOGMatches++
			}

			if showReads {
				fmt.Println(pileup.FormatRecord(rec))
				break
			} else if verbose {
				fmt.Printf("%c%d%c%s depth:%d rank:%d OG:%t\n",
					g.Nts[0][i], rec.Pos+1, read.Nt, silentS, read.Depth,
					rank, matchesOg)
			}
		}
	}
	fmt.Printf("%d/%d are OG matches\n", totalOGMatches, total)
	OR, p := og.IsRemarkable(totalOGMatches, total)
	fmt.Printf("OR=%.2f p=%.5g\n", OR, p)
}

func main() {
	var (
		fasta, orfs string
		recCAS      string
		minDepth    int
		minRatio    float64
		silent      bool
		tc          bool
		reparse     bool
		showReads   bool
		quiet       bool
	)

	flag.StringVar(&fasta, "ref", "", "Reference genome")
	flag.StringVar(&recCAS, "recCA", "", "RecCA genome")
	flag.StringVar(&orfs, "orfs", "", "Reference ORFs")
	flag.IntVar(&minDepth, "min-depth", 4, "Minimum depth")
	flag.Float64Var(&minRatio, "min-ratio", 0, "Minimum ratio of QS to majority")
	flag.BoolVar(&silent, "silent", false, "Require silent")
	flag.BoolVar(&reparse, "reparse", false, "Parse our own .txt.gz format")
	flag.BoolVar(&tc, "tc", false, "Only look at TC")
	flag.BoolVar(&showReads, "show-reads",
		false, "Just show the reads in our pileup format")
	flag.BoolVar(&quiet, "q", false, "Just output OR/p")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)

	if recCAS == "" {
		recCAS = path.Join(os.Getenv("GOPATH"),
			"src/genomics/fasta/recCA.fasta")
	}
	recCA := genomes.LoadGenomes(recCAS, orfs, false)

	var og Outgroup
	og.Init(recCA, g)

	for _, arg := range flag.Args() {
		var pu *pileup.Pileup
		var err error

		if reparse {
			pu, err = pileup.Parse2(arg)
		} else {
			pu, err = pileup.Parse(arg)
		}
		if err != nil {
			log.Fatal(err)
		}
		Compare(pu, g, minDepth, minRatio, silent, tc, showReads, !quiet, &og)
	}
}
