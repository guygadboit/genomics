package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"io"
	"log"
	"os"
	"strings"
)

type Source struct {
	name    string
	index   string
	ntFreq  map[string]float64
	dinFreq map[string]float64
}

func loadFrequency(source Source, numNts int) map[string]float64 {
	ret := make(map[string]float64)

	f := utils.NewFileReader(fmt.Sprintf(
		"../subsequences/output/%s-%d.txt", source.name, numNts))
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

		nt := strings.TrimRight(fields[0], ":")
		ret[nt] = float64(utils.Atoi(fields[1]))
		total += ret[nt]
	}

	for k, _ := range ret {
		ret[k] /= total
	}

	return ret
}

func getSources() []Source {
	root := "/fs/f/genomes/"

	sources := []Source{
		{"Human", root + "human/index", nil, nil},
		{"Cod", root + "cod/index", nil, nil},
		{"DR", root + "bacteria/GCRich/dr_index", nil, nil},
		{"Legionella", root + "bacteria/Legionella/index", nil, nil},
		{"Salmonella", root + "bacteria/Salmonella/index", nil, nil},
		{"Ricksettia", root + "bacteria/Ricksettia/index", nil, nil},
		{"HI", root + "bacteria/ATRich/hi_index", nil, nil},
		{"PA", root + "bacteria/PseudomonasAeruginosa/index", nil, nil},
		{"Listeria", root + "bacteria/Listeria/index", nil, nil},
		{"Streptomyces", root + "bacteria/Streptomyces/index", nil, nil},
		{"Delftia", root + "bacteria/delftia/index", nil, nil},
	}

	for i := 0; i < len(sources); i++ {
		sources[i].ntFreq = loadFrequency(sources[i], 1)
		sources[i].dinFreq = loadFrequency(sources[i], 2)
	}

	return sources
}

/*
	What is the expected frequency of pat based on the frequencies of either
	the individual nts or the dinucleotides of which it is composed?
*/
func expectedFrequency(pat []byte,
	freq map[string]float64, numNts int) float64 {
	ret := 1.0
	for i := 0; i < len(pat)-numNts; i++ {
		ret *= freq[string(pat[i:i+numNts])]
	}
	return ret
}

func count(ins *Insertion, source *Source) (int, float64, float64) {
	var (
		search genomes.BidiIndexSearch
		count  int
	)

	for search.Init(source.index, ins.nts); !search.End(); search.Next() {
		count++
	}

	freq := float64(count) / float64(search.GenomeLength())
	or := freq / expectedFrequency(ins.nts, source.ntFreq, 1)
	or2 := (freq * freq) / expectedFrequency(ins.nts, source.dinFreq, 2)

	return count, or, or2
}

func countInGenomes(insertions []Insertion,
	filters []filterFunc, verbose bool) {
	seen := make(map[string]bool)
	sources := getSources()
	var numProcessed int

	fd, err := os.Create("or-results.txt")
	if err != nil {
		log.Fatal("Can't create results file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	fmt.Fprintf(w, "pattern ")
	for i := 0; i < len(sources); i++ {
		fmt.Fprintf(w, "%s ", sources[i].name)
	}
	fmt.Fprintf(w, "\n")

	filterInsertions(insertions, filters, func(ins *Insertion) {
		// Unique them as we go along
		nts := string(ins.nts)
		there := seen[nts]
		if there {
			return
		}
		seen[nts] = true

		fmt.Fprintf(w, "%s ", string(ins.nts))
		for i := 0; i < len(sources); i++ {
			num, or, or2 := count(ins, &sources[i])
			fmt.Fprintf(w, "%d,%.2f,%.2f ", num, or, or2)
		}
		fmt.Fprintf(w, "\n")

		numProcessed++
		if numProcessed%10 == 0 {
			fmt.Printf("Processed %d\n", numProcessed)
		}

	}, false)
	w.Flush()
}
