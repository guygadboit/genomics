package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/mutations"
	"io"
	"math/rand"
)

type TamperTrial struct {
	runFunc func(genome *genomes.Genomes,
		numMuts int, results chan interface{})
}

func (t *TamperTrial) Run(genome *genomes.Genomes,
	numMuts int, results chan interface{}) {
	t.runFunc(genome, numMuts, results)
}

func (t *TamperTrial) WriteHeadings(w io.Writer) {
	fmt.Fprintln(w, "# Results from a Tamper Trial")
	fmt.Fprintln(w, "name tampered muts_in_sites total_sites total_singles acc")
}

type TamperTrialResult struct {
	SilentInSites
	name             string
	tampered         bool
	detectedTampered bool
}

func (r *TamperTrialResult) Write(w io.Writer) {
	fmt.Fprintln(w, r.name, r.tampered,
		r.totalMuts, r.totalSites, r.totalSingleSites, r.detectedTampered)
}

func TamperTrials(genome *genomes.Genomes, nd *mutations.NucDistro,
	numTrials int, numMuts int, numEdits int, results chan interface{}) {

	classifier := GetClassifier()

	reportProgress := func(n int) {
		fmt.Printf("%s (%d muts) %d/%d trials\n",
			genome.Names[0], numMuts, n, numTrials)
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Filter(0, 0)
		mutant.DeepCopy(0)
		mutations.MutateSilent(mutant, nd, numMuts, 1)

		tampered := rand.Intn(2) == 1
		if tampered {
			Tamper(mutant, RE_SITES, numEdits, numEdits)
		}

		var result TamperTrialResult
		mutant.Combine(genome)
		result.SilentInSites = CountSilentInSites(mutant, RE_SITES, true)
		result.name = genome.Names[0]
		result.tampered = tampered
		result.detectedTampered = classifier.IsTampered(mutant)

		results <- &result

		if i%100 == 0 {
			reportProgress(i)
		}
	}
}
