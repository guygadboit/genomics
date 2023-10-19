/*
	The idea here is for each insertion that you find in some bacteria, output
	the bits either side into a fastq file. We will then bowtie2 them against
	SC2. Any that align might imply that the source of them is contamination--
	but the original assembly included them because they happened to match
	either side of the insertion.
*/
package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"log"
	"os"
)

func outputFastq(insertions []Insertion, filters []filterFunc,
	sources []Source, outName string) {
	fd, err := os.Create(outName)
	if err != nil {
		log.Fatalf("Can't create %s", outName)
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	seen := make(map[string]bool)
	// Cache genomes as we load them here
	cache := make(map[string]*genomes.Genomes)

	filterInsertions(insertions, filters, func(ins *Insertion) {
		// Unique them as we go along
		nts := string(ins.nts)
		there := seen[nts]
		if there {
			return
		}
		seen[nts] = true

		for i := 0; i < len(sources); i++ {
			source := &sources[i]
			var search genomes.BidiIndexSearch

			for search.Init(source.index,
				ins.nts); !search.End(); search.Next() {
				pos, _ := search.Get()

				genome, there := cache[source.fasta]
				if !there {
					genome = genomes.LoadGenomes(source.fasta, "", true)
					cache[source.fasta] = genome
				}

				pre := genome.Slice(0, pos-20, pos)
				after := pos + len(ins.nts)
				post := genome.Slice(0, pos, after+20)
				len := len(pre) + len(post)

				fmt.Fprintf(w, "@%s-%d.%d length=%d\n",
					source.name, ins.id, pos, len)
				fmt.Fprintf(w, "%s%s\n+\n", pre, post)

				for j := 0; j < len; j++ {
					fmt.Fprintf(w, "?")
				}
				fmt.Fprintf(w, "\n")
			}
		}
	}, false)

	w.Flush()
}
