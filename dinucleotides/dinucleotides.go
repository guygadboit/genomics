package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"log"
	"os"
)

type DinucMap map[string]int

// Reject anything with funky ambiguous nts or Ns etc.
func isValid(s string) bool {
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case 'G':
			fallthrough
		case 'A':
			fallthrough
		case 'T':
			fallthrough
		case 'C':
			continue
		default:
			return false
		}
	}
	return true
}

/*
	Count the dinucleotides and return a map of their counts
*/
func FindDinucs(genome *genomes.Genomes, which int) DinucMap {
	ret := make(DinucMap)
	nts := genome.Nts[which]

	for i := 0; i < genome.Length()-1; i++ {
		dn := string(nts[i : i+2])
		if !isValid(dn) {
			continue
		}
		count, _ := ret[dn]
		ret[dn] = count + 1
	}
	return ret
}

func (d DinucMap) Save(fname string) {
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)
	for k, v := range d {
		fmt.Fprintf(w, "%s: %d\n", k, v)
	}
	w.Flush()
}

type Source struct {
	name  string
	fname string
}

func getSources() []Source {
	root := "/fs/f/genomes/"
	return []Source{
		{"Viruses", root + "viruses/mega/mega.fasta"},
		{"Bat", root + "bat/myotis_davidii/" +
			"GCF_000327345.1_ASM32734v1_genomic.fna.gz"},
		{"Human", root + "human/GRCh38_latest_genomic.fna.gz"},
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
		{"Insertions", "../fasta/InsertionsNotFromWH1.fasta"},
	}
}

func main() {
	sources := getSources()

	for i := 0; i < len(sources); i++ {
		source := sources[i]
		g := genomes.LoadGenomes(source.fname, "", true)
		dm := FindDinucs(g, 0)

		outName := fmt.Sprintf("%s-dinucs.txt", source.name)
		dm.Save(outName)
		fmt.Printf("Wrote %s\n", outName)
	}
}