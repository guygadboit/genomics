package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
)

func testMutations(genome *genomes.Genomes) {
	fmt.Printf("Loaded %d genomes length %d\n",
		genome.NumGenomes(), genome.Length())

	genome.Save("B52", "B52-test.fasta", 0)
	fmt.Printf("Saved as B52-test.fasta\n")

	nd := mutations.NewNucDistro(mutations.NewGenomeIterator(genome),
		mutations.NT_ALPHABET)

	var mutant *genomes.Genomes
	for {
		mutant = genome.Clone()
		mutations.MutateSilent(mutant, nd, 1, 700)
		count, maxLength, unique, interleaved, _, _ :=
			FindRestrictionMap(mutant)
		if unique && maxLength < 8000 {
			fmt.Println(count, maxLength, unique, interleaved)
			break
		}
	}
	mutant.Save("Mutant", "B52-mutated.fasta", 0)
	fmt.Printf("Saved as B52-mutated.fasta\n")
}

func testTamper(genome *genomes.Genomes) {
	num := Tamper(genome, RE_SITES, 10, 10)
	fmt.Printf("Tampered with %d sites\n", num)

	genome.Save("Mutant", "B52-mutated.fasta", 0)
	fmt.Printf("Saved as B52-mutated.fasta\n")
}

func testAlternatives(genome *genomes.Genomes) {
	var env genomes.Environment
	env.Init(genome, 300, 6, 0)

	alternatives := env.FindAlternatives(6, true)
	fmt.Println(alternatives)
}

func testCachedSearch(genome *genomes.Genomes) {
	var cs CachedSearch

	for i := 0; i < 3; i++ {
		fmt.Printf("Starting search\n")
		cs.Init(genome, RE_SITES)
		for {
			pos, site := cs.Iter()
			if cs.End() {
				break
			}

			fmt.Printf("%s at %d\n", string(site.pattern), pos)
		}
	}
}

func testTranslate(genome *genomes.Genomes) {
	var it genomes.CodonIter
	it.Init(genome, 0)

	for {
		pos, codon, aa, err := it.Next()
		if err != nil {
			break
		}
		fmt.Printf("%d: %s %c\n", pos, codon, aa)
	}
	fmt.Printf("\n")
}

func Test() {
	genome := genomes.LoadGenomes("../fasta/B52.fasta",
		"../fasta/B52.orfs", false)
	// testCachedSearch(genome)
	// testMutations(genome)
	// testAlternatives(genome)
	// testTamper(genome)
	testTranslate(genome)
}
