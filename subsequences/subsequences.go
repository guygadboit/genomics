package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"log"
	"os"
	"sort"
	"sync"
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

type Source struct {
	name  string
	fname string
}

func getSources() []Source {
	root := "/fs/f/genomes/"
	return []Source{
		{"Bat", root + "bat/myotis_davidii/" +
			"GCF_000327345.1_ASM32734v1_genomic.fna.gz"},
		{"Human", root + "human/GRCh38_latest_genomic.fna.gz"},
		{"RaccoonDog", root + "raccoon_dog/" +
			"GCF_905146905.1_NYPRO_anot_genome_genomic.fna.gz"},
		{"Pangolin", root + "pangolin/" +
			"GCF_014570535.1_YNU_ManJav_2.0_genomic.fna.gz"},
		{"Viruses", root + "viruses/mega/mega.fasta"},
	}
}

func writeResults(source Source, ss SubSequences, length int) {
	fname := fmt.Sprintf("%s-%d.txt", source.name, length)
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create results file")
	}
	defer fd.Close()
	w := bufio.NewWriter(fd)

	for i := 0; i < len(ss); i++ {
		fmt.Fprintf(w, "%s: %d\n", string(ss[i].nts), ss[i].count)
	}
	w.Flush()

	fmt.Printf("Wrote %s\n", fname)
}

func findAll(length int) {
	sources := getSources()
	var wg sync.WaitGroup
	maxThreads := 2
	pending := 0

processing:
	for {
		for i := 0; i < maxThreads; i++ {
			if pending < len(sources) {
				wg.Add(1)
				source := sources[pending]
				pending++
				go func() {
					g := genomes.LoadGenomes(source.fname, "", true)
					ss, _ := FindSubSequences(g, 0, length)
					writeResults(source, ss, length)
					wg.Done()
				}()
			} else {
				wg.Wait()
				break processing
			}
		}
		wg.Wait()
	}
}

func main() {
	findAll(6)
}
