package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"log"
	"os"
	"sort"
)

type SubSequence struct {
	nts   string
	count int
}

type SubSequences []SubSequence

func (s SubSequences) Len() int {
	return len(s)
}

func (s SubSequences) Less(i, j int) bool {
	return s[i].count < s[j].count
}

func (s SubSequences) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

/*
	Find all the subsequences and their lengths. Return a sorted list of each
	one, how many times it occurs, and the total.
*/
func FindSubSequences(genome *genomes.Genomes,
	which int, length int) (SubSequences, int) {
	nts := genome.Nts[which]
	subseqs := make(map[string]int)
	var total int

searching:
	for i := 0; i < len(nts)-length; i++ {
		s := nts[i : i+length]

		// Ignore any with Ns or other crap in them
		for j := 0; j < length; j++ {
			switch s[j] {
			case 'G':
				break
			case 'A':
				break
			case 'T':
				break
			case 'C':
				break
			default:
				continue searching
			}
		}

		ss := string(s)
		count, _ := subseqs[ss]
		subseqs[ss] = count + 1
		total++

		if i%100000000 == 0 {
			done := float64(i*100) / float64(len(nts)-length)
			fmt.Printf("%.2f%% done\n", done)
		}
	}

	ret := make(SubSequences, 0, len(subseqs))
	for k, v := range subseqs {
		ret = append(ret, SubSequence{k, v})
	}

	sort.Sort(ret)
	return ret, total
}

func main() {
	fname := "/fs/f/genomes/bat/myotis_davidii/" +
		"GCA_000327345.1_ASM32734v1_genomic.fna.gz"
	//fname := "/fs/f/genomes/human/GRCh38_latest_genomic.fna.gz"
	fmt.Println("Loading genome...")
	g := genomes.LoadGenomes(fname, "", true)

	fmt.Printf("%d nts. Finding subsequences...\n", len(g.Nts[0]))
	ss, total := FindSubSequences(g, 0, 6)

	fmt.Printf("Total: %d\n", total)

	fd, err := os.Create("results.txt")
	if err != nil {
		log.Fatal("Can't create results file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	for i := 0; i < len(ss); i++ {
		count := ss[i].count
		rate := float64(count*1000) / float64(total)
		fmt.Fprintf(w, "%s: %d (%.4f permille)\n",
			string(ss[i].nts), count, rate)
	}

	fmt.Printf("Wrote results.txt\n")
}
