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
	name   string
	index  string
	ntFreq map[byte]float64
}

func loadFrequency(source Source) map[byte]float64 {
	ret := make(map[byte]float64)

	f := utils.NewFileReader(fmt.Sprintf(
		"../subsequences/output/%s-1.txt", source.name))
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

func getSources() []Source {
	root := "/fs/f/genomes/"

	sources := []Source{
		{"Human", root + "human/index", nil},
		{"DR", root + "bacteria/GCRich/dr_index", nil},
		{"Legionella", root + "bacteria/Legionella/index", nil},
		{"Salmonella", root + "bacteria/Salmonella/index", nil},
		{"Ricksettia", root + "bacteria/Ricksettia/index", nil},
		{"HI", root + "bacteria/ATRich/hi_index", nil},
		{"PA", root + "bacteria/PseudomonasAeruginosa/index", nil},
		{"Listeria", root + "bacteria/Listeria/index", nil},
	}

	for i := 0; i < len(sources); i++ {
		sources[i].ntFreq = loadFrequency(sources[i])
	}

	return sources
}

func expectedFrequency(pat []byte, ntFreq map[byte]float64) float64 {
	ret := 1.0
	for i := 0; i < len(pat); i++ {
		ret *= ntFreq[pat[i]]
	}
	return ret
}

func count(ins *Insertion, source *Source) (int, float64) {
	var search genomes.IndexSearch

	var count int

	for search.Init(source.index, ins.nts); !search.End(); search.Next() {
		count++
	}

	freq := float64(count) / float64(search.GenomeLength())
	or := freq / expectedFrequency(ins.nts, source.ntFreq)

	return count, or
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
			num, or := count(ins, &sources[i])
			fmt.Fprintf(w, "%d,%.2f ", num, or)
		}
		fmt.Fprintf(w, "\n")

		numProcessed++
		if numProcessed%10 == 0 {
			fmt.Printf("Processed %d\n", numProcessed)
		}

	}, false)
	w.Flush()
}
