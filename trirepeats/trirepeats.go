package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"log"
	"os"
)

type Source struct {
	name  string
	fname string
}

func getSources() []Source {
	root := "/fs/f/genomes/"

	return []Source{
		{"PA", root + "bacteria/PseudomonasAeruginosa/" +
			"PseudomonasAeruginosaComplete.fasta"},
		{"CC", root + "bacteria/GCRich/CaulobacterCrescentus.fasta"},
		{"DR", root + "bacteria/GCRich/DeinococcusRadiodurans.fasta"},
		{"Salmonella", root + "bacteria/Salmonella/Salmonella.fasta"},
		{"Listeria", root + "bacteria/Listeria/ListeriaInnocua.fasta"},
		{"Ricksettia", root + "bacteria/Ricksettia/Ricksettia.fasta"},
		{"Legionella", root + "bacteria/Legionella/Legionella.fasta"},
		{"Cod", root + "cod/cod.fasta.gz"},
		{"Caulobacter", root + "bacteria/GCRich/CaulobacterCrescentus.fasta"},
		{"Deinococcus", root + "bacteria/GCRich/DeinococcusRadiodurans.fasta"},
		{"Haemophilus", root + "bacteria/ATRich/HaemophilusInfluenzae.fasta"},
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
		{"HCoVs", "../fasta/HCoVs.fasta"},
		{"Insertions", "../fasta/CombinedInsertions.fasta"},
	}
}

type TriRepeatMap struct {
	genome *genomes.Genomes
	repeats map[string]int
}

/*
	Return all 64 possible trinucleotides in a sorted array
*/
func trinucleotides() []string {
	ret := make([]string, 0, 64)
	nts := []byte{'A', 'C', 'G', 'T'}
	pat := make([]byte, 3)

	for i := 0; i < len(nts); i++ {
		for j := 0; j < len(nts); j++ {
			for k := 0; k < len(nts); k++ {
				pat[0] = nts[i]
				pat[1] = nts[j]
				pat[2] = nts[k]
				ret = append(ret, string(pat))
			}
		}
	}
	return ret
}

func (m *TriRepeatMap) Init(genome *genomes.Genomes) {
	m.genome = genome
	keys := trinucleotides()

	m.repeats = make(map[string]int)
	for i := 0; i < len(keys); i++ {
		m.repeats[keys[i]] = 0
	}
}

func (m *TriRepeatMap) Save(fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	total := float64(m.genome.Length())

	keys := trinucleotides()
	for i := 0; i < len(keys); i++ {
		k := keys[i]
		freq := float64(m.repeats[k] * 1000) / total
		fmt.Fprintf(w, "%s: %f\n", k, freq)
	}

	w.Flush()
}

func CountTriRepeats(genome *genomes.Genomes) TriRepeatMap {
	var ret TriRepeatMap
	ret.Init(genome)

	for i := 0; i < genome.Length()-6; i++ {
		pat := string(genome.Nts[0][i : i+6])
		if pat[:3] != pat[3:] {
			continue
		}
		triNt := pat[:3]
		count := ret.repeats[triNt]
		ret.repeats[triNt] = count + 1
	}
	return ret
}

func main() {
	sources := getSources()
	for i := 0; i < len(sources); i++ {
		g := genomes.LoadGenomes(sources[i].fname, "", true)
		tm := CountTriRepeats(g)
		fname := fmt.Sprintf("%s-3nt-repeats.txt", sources[i].name)
		tm.Save(fname)
		fmt.Printf("Wrote %s\n", fname)
	}
}
