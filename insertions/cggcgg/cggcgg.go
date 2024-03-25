package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"regexp"
	"os"
)

func search(needle, haystack []byte) int {
	var count int
outer:
	for i := 0; i < len(haystack)-len(needle)+1; i++ {
		for j := 0; j < len(needle); j++ {
			if needle[j] != haystack[i+j] {
				continue outer
			}
		}
		count++
	}
	return count
}

func searchAligned(needle, haystack []byte) int {
	var count int
outer:
	for i := 0; i < len(haystack)-len(needle)+1; i += 3 {
		for j := 0; j < len(needle); j++ {
			if needle[j] != haystack[i+j] {
				continue outer
			}
		}
		count++
	}
	return count
}

func countRR(g *genomes.Genomes, reverse bool, w *bufio.Writer) (int, int) {
	var cggTotal, otherRRTotal int
	needle := []byte("CGGCGG")
	if reverse {
		needle = utils.ReverseComplement(needle)
	}
	re := regexp.MustCompile("_[^_]+$")

	for i := 0; i < g.NumGenomes(); i++ {
		for offset := 0; offset < 3; offset++ {
			nts := g.Nts[i]
			if reverse {
				nts = utils.ReverseComplement(nts)
			}
			nts = nts[offset:]
			cggCount := searchAligned(needle, nts)
			cggTotal += cggCount

			if cggCount != 0 {
				fmt.Printf("%s: found %s %d times in %s\n",
					g.Names[i], string(needle), cggCount, string(nts))
			}

			protein := genomes.TranslateAlignedShort(nts)
			otherRR := search([]byte("RR"), protein) - cggCount
			if otherRR < 0 {
				otherRR = 0
			}

			if otherRR != 0 {
				fmt.Printf("%s: found non-%s RR %d times in %s (%d: %s)\n",
					g.Names[i], string(needle),
					otherRR, string(nts), offset, string(protein))
				otherRRTotal += otherRR
			}

			var direction string
			if reverse {
				direction = "3'->5'"
			} else {
				direction = "5'->3'"
			}

			shortName := string(re.ReplaceAll([]byte(g.Names[i]), []byte{}))
			url := fmt.Sprintf("https://cov-spectrum.org/explore/"+
				"World/AllSamples/AllTimes/variants?nucInsertions=%s%%3A%s",
				shortName, string(g.Nts[i]))

			if cggCount != 0 || otherRR != 0 {
				fmt.Fprintf(w, "%s:%s,%s,%d,%d,%d,%s,%s\n", shortName,
					string(g.Nts[i]), direction,
					cggCount, otherRR, offset, string(protein), url)
			}
		}
	}
	return cggTotal, otherRRTotal
}

func main() {
	g := genomes.LoadGenomes("../../fasta/SplitInsertions.fasta", "", false)
	f, err := os.Create("results.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	fmt.Fprintln(w, "Insertion,Direction,CGGCGG Count,"+
		"Other RR Count,Offset,Translation,URL")

	var yes, no int
	a, b := countRR(g, false, w)
	yes += a
	no += b

	a, b = countRR(g, true, w)
	yes += a
	no += b

	total := yes + no
	fmt.Printf("%d %d %.2f%%\n", yes, yes+no, float64(yes*100)/float64(total))

	w.Flush()
	fmt.Printf("Wrote %s\n", "results.csv")
}
