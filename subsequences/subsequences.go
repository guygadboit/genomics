package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"io"
	"log"
	"os"
	"sort"
	"strings"
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
	// root := "/fs/f/genomes/"

	sources := []Source{
		{"Yeast", "../fasta/Yeast.fasta"},
	}
	/*
			{"mRNAVaccines", "../fasta/mRNAVaccines.fasta.gz"},
			{"Viruses", root + "viruses/mega/mega.fasta"},
			{"Bat", root + "bat/myotis_davidii/" +
				"GCF_000327345.1_ASM32734v1_genomic.fna.gz"},
			{"RaccoonDog", root + "raccoon_dog/" +
				"GCF_905146905.1_NYPRO_anot_genome_genomic.fna.gz"},
			{"Pangolin", root + "pangolin/" +
				"GCF_014570535.1_YNU_ManJav_2.0_genomic.fna.gz"},
			{"Rabbit", root + "rabbit/" +
				"GCF_009806435.1_UM_NZW_1.0_genomic.fna.gz"},
			{"Streptomyces", root + "bacteria/Streptomyces/" +
				"GCF_000009765.2_ASM976v2_genomic.fna.gz"},
			{"Pig", root + "pig/" +
				"GCF_000003025.6_Sscrofa11.1_genomic.fna.gz"},
			{"Mouse", root + "mouse/" +
				"GCF_000001635.27_GRCm39_genomic.fna.gz"},
			{"Cod", root + "cod/" + "GCF_902167405.1_gadMor3.0_genomic.fna.gz"},
		}
	*/

	/*
		bacteria := []string{
			"Delftia", "DR", "Legionella",
			"Salmonella", "Ricksettia", "HI",
			"PA", "Listeria", "Streptomyces",
			"StrepPyogenes", "StrepPneum", "Mycoplasma",
			"Brucella", "OT", "RP",
		}

		for i := 0; i < len(bacteria); i++ {
			name := bacteria[i]
			sources = append(sources, Source{name,
				fmt.Sprintf("../fasta/bacteria/%s/%s.fasta.gz", name, name)})
		}
	*/

	return sources
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

func findAll(sources []Source, length int) {
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

func findPattern(sources []Source, pattern []byte) {
	var s genomes.Search
	// fcs := []byte("CTCCTCGGCGGG")

	// OR of the below occurring in Yeast based on individual nt frequencies is
	// 2402.6004106240753 with a p of 1.0 (binomial test)
	// fcs := []byte("TTCTCCTCGGCGGGCA")

	for i := 0; i < len(sources); i++ {
		var count int
		source := sources[i]
		g := genomes.LoadGenomes(source.fname, "", true)
		for s.Init(g, 0, pattern, 0.0); !s.End(); s.Next() {
			count++
		}
		total := g.Length()
		freq := float64(count) / float64(total)
		fmt.Printf("%s: %d/%d (%g)\n", source.name, count, total, freq)

		expected := loadExpected(source)

		expFreq := expectedFrequency(pattern, expected)

		fmt.Printf("Expected frequency: %g OR %f\n", expFreq, freq/expFreq)
		fmt.Printf("Contingency table: %d %d 1 %d\n",
			count, total, int(1/expFreq))
	}
}

func loadExpected(source Source) map[byte]float64 {
	ret := make(map[byte]float64)

	f := utils.NewFileReader(fmt.Sprintf("output/%s-1.txt", source.name))
	defer f.Close()

	var total float64
loop:
	for {
		line, err := f.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)
		fields := strings.Fields(line)

		nt := []byte(strings.TrimRight(fields[0], ":"))[0]
		ret[nt] = float64(utils.Atoi(fields[1]))
		total += ret[nt]
	}

	for k, _ := range ret {
		ret[k] /= total
	}

	return ret
}

func expectedFrequency(pat []byte, ntFreq map[byte]float64) float64 {
	ret := 1.0
	for i := 0; i < len(pat); i++ {
		ret *= ntFreq[pat[i]]
	}
	return ret
}

func montecarlo(length int, nTrials int, index string, verbose bool) int {
	sources := getSources()
	nFound := 0

	testGenomes := make([]*genomes.Genomes, len(sources))
	expected := make([]map[byte]float64, len(sources))

	for i, source := range sources {
		testGenomes[i] = genomes.LoadGenomes(source.fname, "", true)
		expected[i] = loadExpected(source)
	}

	initSearch := func(haystack *genomes.Genomes,
		pat []byte) genomes.SearchIf {
		if index != "" {
			var is genomes.IndexSearch
			is.Init(index, pat)
			return &is
		} else {
			var rs genomes.Search
			rs.Init(haystack, 0, pat, 0.0)
			return &rs
		}
	}

	for i := 0; i < nTrials; i++ {
		pat := utils.RandomNts(length)
		var count int
		for j := 0; j < len(testGenomes); j++ {
			g := testGenomes[j]
			for s := initSearch(g, pat); !s.End(); s.Next() {
				count++
			}
			freq := float64(count) / float64(g.Length())
			or := freq / expectedFrequency(pat, expected[j])
			if count > 0 {
				nFound++
				if verbose {
					fmt.Printf("%d/%d %s: %.4f %s\n", i, nTrials,
					sources[j].name, or, string(pat))
				}
			}
		}
	}
	fmt.Printf("Random sequences of length %d occurred %d/%d (%g%%)\n",
		length, nFound, nTrials, float64(nFound*100)/float64(nTrials))

	return nFound
}

func main() {
	var pat string
	var length int
	var source string
	var monte int
	var index string
	var verbose bool

	// Note: in order to look for a pat first you need to run with length 1 and
	// copy the output into the output directory (where you might want to add
	// it to git)
	flag.StringVar(&pat, "p", "", "Pattern of interest")
	flag.StringVar(&source, "s", "", "Source fasta file")
	flag.StringVar(&index, "index", "", "Directory of index (optional)")
	flag.IntVar(&length, "l", 0, "Length of subsequences to find freq of")
	flag.IntVar(&monte, "m", 0, "Montecarlo iterations")
	flag.BoolVar(&verbose, "v", false, "Verbose")

	flag.Parse()

	var sources []Source
	if source == "" {
		sources = getSources()
	} else {
		sources = make([]Source, 1)
		name := utils.BasicName(source)
		sources[0] = Source{name, source}
	}

	if pat != "" {
		findPattern(sources, []byte(pat))
		return
	}

	if monte != 0 {
		montecarlo(length, monte, index, verbose)
		return
	}

	if length != 0 {
		findAll(sources, length)
	}
}
