/*
Run this on a pileup generated with e.g.:

$ samtools mpileup sorted-output.sam

And the reference genome that the sam file was originally generated from. It
outputs an alignment of the two in fasta format
*/
package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/pileup"
	"genomics/utils"
	"log"
	"math"
	"os"
	"strings"
)

func ConsensusSubsequence(p *pileup.Pileup, start, end int) string {
	ret := ""
	for pos := start; pos < end; pos++ {
		record := p.Get(pos)
		if record == nil {
			ret += "N"
			continue
		}
		ret += fmt.Sprintf("%c", record.Reads[0].Nt)
	}
	return ret
}

func Match(pileup *pileup.Pileup, pattern []byte,
	pos int, minDepth int, tolerance float64) bool {
	allowedFails := int(math.Floor(tolerance) * float64(len(pattern)))

	fails := 0
outer:
	for i, c := range pattern {
		record := pileup.Get(pos + i)
		if record == nil {
			fails++
			continue
		}
		for _, read := range record.Reads {
			if read.Nt == c && read.Depth >= minDepth {
				continue outer
			}
			fails++
		}
		if fails > allowedFails {
			break
		}
	}
	return fails < allowedFails
}

func ParseMatchExpr(matchExpr string) (int, []byte) {
	fields := strings.Split(matchExpr, ":")
	return utils.Atoi(fields[0]) - 1, []byte(fields[1])
}

func main() {
	var (
		reference, output    string
		verbose, veryVerbose bool
		justShow             bool
		showPosS             string
		subseq               string
		matchExpr            string
		matchTol             float64
		matchMinDepth        int
		reparse              bool
	)

	flag.StringVar(&reference, "ref", "", "Reference genome")
	flag.StringVar(&output, "o", "output.fasta", "Output name")
	flag.BoolVar(&verbose, "v", false, "verbose")
	flag.BoolVar(&veryVerbose, "vv", false, "very verbose")
	flag.BoolVar(&justShow, "show", false, "Just show counts in the pileup")
	flag.StringVar(&showPosS, "pos", "", "Positions to show counts for")
	flag.StringVar(&subseq, "consensus-subseq",
		"", "Show consensus for a subsequence")
	flag.StringVar(&matchExpr, "match",
		"", "Match an expression of the form pos:pattern above min depth")
	flag.Float64Var(&matchTol, "match-tol", 0.2, "Match tolerance")
	flag.IntVar(&matchMinDepth, "match-min-depth", 6, "Match min depth")
	flag.BoolVar(&reparse, "reparse", false, "Parse our own previous show"+
		"output rather than an mpileup file")
	flag.Parse()

	if len(flag.Args()) != 1 {
		flag.PrintDefaults()
		return
	}

	var pu *pileup.Pileup
	var err error
	if reparse {
		pu, err = pileup.Parse2(flag.Args()[0])
	} else {
		pu, err = pileup.Parse(flag.Args()[0])
	}
	if err != nil {
		log.Fatal(err)
	}

	var showPositions []int
	if showPosS != "" {
		justShow = true
		showPositions = utils.ParseInts(showPosS, ",")
		for i, v := range showPositions {
			showPositions[i] = v - 1
		}
	}

	if justShow {
		pu.Show(showPositions)
		return
	}

	if subseq != "" {
		ss := utils.ParseInts(subseq, ":")
		fmt.Println(ConsensusSubsequence(pu, ss[0]-1, ss[1]))
		return
	}

	if matchExpr != "" {
		pos, pattern := ParseMatchExpr(matchExpr)
		matched := Match(pu, pattern, pos, matchMinDepth, matchTol)
		cs := ConsensusSubsequence(pu, pos, pos+len(pattern))
		var matchS string
		if matched {
			matchS = "matched"
		}
		fmt.Println(cs, matchS)
		if !matched {
			os.Exit(-1)
		} else {
			return
		}
	}

	g := genomes.LoadGenomes(reference, "", false)

	// The second genome will be what we find in the pileup
	g = g.Filter(0, 0)
	g.DeepCopy(1)
	g.Names[1] = "Pileup"

	nts := g.Nts[1]
	var j int
	for i := 0; i < g.Length(); i++ {
		if j < len(pu.Records) && pu.Records[j].Pos == i {
			nts[i] = pu.Records[j].Reads[0].Nt

			doPrint := veryVerbose || (verbose && nts[i] != g.Nts[0][i])

			if doPrint {
				for _, r := range pu.Records[j].Reads {
					fmt.Printf("%d%c %d\n", i+1, r.Nt, r.Depth)
				}
			}

			j++
		} else {
			nts[i] = '-'
		}
	}

	if output != "" {
		g.SaveMulti(output)
		fmt.Printf("Wrote %s\n", output)
	}
}
