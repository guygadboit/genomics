package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sync"
	"genomics/genomes"
	"genomics/mutations"
)

type TrialResult interface {
	Write(w io.Writer)
}

type Trial interface {
	WriteHeadings(w io.Writer)
	Run(genome *genomes.Genomes, numMuts int, results chan interface{})
}

func loadGenomes(fnames []string) []*genomes.Genomes {
	g := make([]*genomes.Genomes, len(fnames))
	for i := 0; i < len(fnames); i++ {
		g[i] = genomes.LoadGenomes(
			"../fasta/"+fnames[i]+".fasta",
			"../fasta/"+fnames[i]+".orfs",
			false)
	}
	return g
}

func findNucDistro(g []*genomes.Genomes) *mutations.NucDistro {
	nd := mutations.NewNucDistro(nil)
	for i := 0; i < len(g); i++ {
		nd.Count(g[i])
	}
	return nd
}

/*
Return an array of ints for how many muts to apply, which either numMuts
for everything, or the number of silent muts there are between each genome
and WH1.
*/
func findMutsPerGenome(fnames []string, numMuts int) []int {
	mutsPerGenome := make([]int, len(fnames))

	for i := 0; i < len(fnames); i++ {
		if numMuts != 0 {
			mutsPerGenome[i] = numMuts
		} else {
			g := genomes.LoadGenomes(fmt.Sprintf("WH1-%s.fasta",
				fnames[i]), "../fasta/WH1.orfs", false)
			mutsPerGenome[i], _ = mutations.CountMutations(g)
		}
	}

	return mutsPerGenome
}

func writeParams(w io.Writer, nTrials, nMuts, nEdits int) {
	fmt.Fprintf(w, "# Trials: %d Muts: %d (0 means auto) Edits: %d\n",
		nTrials, nMuts, nEdits)
}

func main() {
	var nTrials, nMuts, nThreads, nEdits int
	var test, countSites bool
	var trialType string

	flag.IntVar(&nTrials, "n", 10000, "Number of trials")
	flag.IntVar(&nMuts, "m", 0, "Number of mutations (0 means auto)")
	flag.IntVar(&nThreads, "p", 1, "Number of threads")
	flag.BoolVar(&test, "t", false, "Just do some self-tests")
	flag.BoolVar(&countSites, "c", false, "Count mutations per site etc.")
	flag.StringVar(&trialType, "trial", "spacing", "Which trials to run")
	flag.IntVar(&nEdits, "edits", 3, "Number of sites to move")
	flag.Parse()

	if test {
		Test()
		return
	}

	var c Classifier
	c.Init()
	c.Trial(100)
	return

	fnames := []string{
		"RpYN06",
		"BtSY2",
		"ChimericAncestor",
		"BANAL-20-236",
		"BANAL-20-52",
		"BANAL-20-103",
		"RaTG13",
	}

	g := loadGenomes(fnames)
	nd := findNucDistro(g)
	nd.Show()

	// How many silent muts to apply per genome? If they set 0 that means
	// "auto" so use the same number as there are between that genome and WH1.
	mutsPerGenome := findMutsPerGenome(fnames, nMuts)

	// Construct the trial objects
	spacingTrial := SpacingTrial{
		func(genome *genomes.Genomes, numMuts int, results chan interface{}) {
			SpacingTrials(genome, nd, nTrials/nThreads,
				numMuts, countSites, results)
		}}

	tamperTrial := TamperTrial{
		func(genome *genomes.Genomes, numMuts int, results chan interface{}) {
			TamperTrials(genome, nd,
				nTrials/nThreads, numMuts, nEdits, results)
		}}

	trials := map[string]Trial{
		"spacing": &spacingTrial,
		"tamper":  &tamperTrial,
	}

	trial := trials[trialType]

	fd, err := os.Create("results.txt")
	if err != nil {
		log.Fatal("Can't create results file")
	}
	defer fd.Close()

	resultsWriter := bufio.NewWriter(fd)
	writeParams(resultsWriter, nTrials, nMuts, nEdits)

	trial.WriteHeadings(resultsWriter)
	results := make(chan interface{}, 1000)

	if trialType == "tamper" {
		// Write the reference values into the results file
		for i := 0; i < len(fnames); i++ {
			CountSilentInSitesReference(fnames[i], RE_SITES, results)
		}
	}

	var wg sync.WaitGroup

	// Cut the work up unto nThreads pieces, all writing their results to a
	// single channel. Each thread will do a portion of the tests but for all
	// genomes
	for i := 0; i < nThreads; i++ {
		wg.Add(len(g))

		go func() {
			for j := 0; j < len(g); j++ {
				trial.Run(g[j], mutsPerGenome[j], results)
				wg.Done()
			}
		}()
	}

	// Keep reading out of the results channel and writing to the results file
	// until everyone has finished (we will write to stop once wg has
	// completed)
	stop := make(chan bool)
	go func(stop chan bool) {
	loop:
		for {
			select {
			case r := <-results:
				trialResult := r.(TrialResult)
				trialResult.Write(resultsWriter)
			case <-stop:
				break loop
			}
		}
		resultsWriter.Flush()
		fmt.Println("Wrote results.txt")
	}(stop)

	wg.Wait()
	stop <- true
}
