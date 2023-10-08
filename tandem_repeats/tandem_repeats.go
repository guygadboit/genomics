package main

import (
	"genomics/genomes"
	"fmt"
)

func findTandemRepeats(genome *genomes.Genomes,
	length int, verbose bool) []int {
	ret := make([]int, 0)
	nts := genome.Nts[0]
searching:
	for i := 0; i < genome.Length() - length*2; i++ {
		for j := 0; j < length; j++ {
			if nts[i+j] != nts[i+length+j] {
				continue searching
			}
		}
		ret = append(ret, i)
		if verbose {
			fmt.Printf("%d (%d): %s\n", i, length, string(nts[i:i+length]))
		}
	}
	return ret
}

func main() {
	root := "/fs/f/genomes/bacteria/"

	fnames := []string{
		"./GCRich/CaulobacterCrescentus.fasta",
		"./GCRich/DeinococcusRadiodurans.fasta",
		"./ATRich/HaemophilusInfluenzae.fasta",
		"./Salmonella/Salmonella.fasta",
		"./Listeria/ListeriaInnocua.fasta",
		"./Ricksettia/Ricksettia.fasta",
		"./Legionella/Legionella.fasta",
		"./PseudomonasAeruginosa/PseudomonasAeruginosaComplete.fasta",
		"./delftia/delftia.fasta.gz",
		"./delftia/delftia.fasta.gz",
		"./Streptomyces/Streptomyces.fasta.gz",
	}

	for i := 0; i < len(fnames); i++ {
		g := genomes.LoadGenomes(root + fnames[i], "", true)
		for length := 12; length <= 60; length++ {
			fmt.Printf("%s\n", fnames[i])
			findTandemRepeats(g, length, true)
		}
	}
}
