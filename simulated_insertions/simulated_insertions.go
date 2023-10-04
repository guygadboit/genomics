package main

import (
	"genomics/genomes"
	"fmt"
	"bufio"
	"os"
	"log"
	"math/rand"
)

func randomInsertions(genome *genomes.Genomes,
	length int, count int, fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create %s\n", fname)
	}
	defer f.Close()

	fp := bufio.NewWriter(f)


	for i := 0; i < count; i++ {
		start := rand.Intn(genome.Length() - length)
		pattern := genome.Nts[0][start:start+length]
		fmt.Fprintf(fp, "%d ins_0000:%s (2 seqs)\n", i, pattern)
	}

	fp.Flush()
}

func main() {
	g := genomes.LoadGenomes(os.Args[1], "", true)
	randomInsertions(g, 14, 100000, os.Args[2])
}
