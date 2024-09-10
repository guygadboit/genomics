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
Count the (di)nucleotides and return a map of their counts. It's
dinucleotides if len is 2.
*/
func FindDinucs(genome *genomes.Genomes, which int, len int) DinucMap {
	ret := make(DinucMap)
	nts := genome.Nts[which]

	for i := 0; i < genome.Length()-len+1; i++ {
		dn := string(nts[i : i+len])
		if !isValid(dn) {
			continue
		}
		ret[dn]++
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
	bacRoot := root + "bacteria/"
	return []Source{

		/*
			{"Insertions", "../fasta/InsertionsNotFromWH1OrHuman.fasta"},
			{"MaybeBac", "../fasta/MaybeBacCombined.fasta"},
			{"PA", "../fasta/bacteria/PA/PA.fasta.gz"},
			{"CC", root + "bacteria/GCRich/CaulobacterCrescentus.fasta.gz"},
			{"DR", root + "bacteria/GCRich/DeinococcusRadiodurans.fasta.gz"},
			{"HI", root + "bacteria/ATRich/HaemophilusInfluenzae.fasta.gz"},
			{"Salmonella", root + "bacteria/Salmonella/Salmonella.fasta.gz"},
			{"Listeria", root + "bacteria/Listeria/ListeriaInnocua.fasta.gz"},
			{"Ricksettia", root + "bacteria/Ricksettia/Ricksettia.fasta.gz"},
			{"Legionella", root + "bacteria/Legionella/Legionella.fasta.gz"},
			{"mRNAVaccines", "../fasta/mRNAVaccines.fasta.gz"},
		*/
		{"Insertions", "../fasta/CombinedInsertions.fasta"},
		{"HI", bacRoot + "ATRich/HaemophilusInfluenzae.fasta.gz"},
		{"AVisc", bacRoot + "AVisc/AVisc.fasta.gz"},
		{"ANaesl", bacRoot + "ANaesl/ANaesl.fasta.gz"},
		{"AIsrael", bacRoot + "AIsrael/AIsrael.fasta.gz"},
		{"Porphyromonas", bacRoot + "Porphyromonas/Porphyromonas.fasta.gz"},
		{"AActinom", bacRoot + "AActinom/AActinom.fasta.gz"},
		{"TForsyth", bacRoot + "TForsyth/TForsyth.fasta.gz"},
		{"Treponema", bacRoot + "Treponema/Treponema.fasta.gz"},
		{"HCoVs", "../fasta/HCoVs.fasta"},
		/*
		{"StrepPyogenes", bacRoot + "StrepPyogenes/StrepPyogenes.fasta.gz"},
		{"StrepPneum", bacRoot + "StrepPneum/StrepPneum.fasta.gz"},
		{"Mycoplasma", bacRoot + "Mycoplasma/Mycoplasma.fasta.gz"},
		{"OT", bacRoot + "OT/OT.fasta.gz"},
		{"Porphyromonas", bacRoot + "Porphyromonas/Porphyromonas.fasta.gz"},
		{"AActinom", bacRoot + "AActinom/AActinom.fasta.gz"},
		{"TForsyth", bacRoot + "TForsyth/TForsyth.fasta.gz"},
		{"Treponema", bacRoot + "Treponema/Treponema.fasta.gz"},
		{"BactFragilis", bacRoot + "BactFragilis/BactFragilis.fasta.gz"},
		*/
		/*
			{"Cod", root + "cod/cod.fasta.gz"},
				{"Caulobacter", root + "bacteria/GCRich/CaulobacterCrescentus.fasta"},
				{"Deinococcus", root + "bacteria/GCRich/DeinococcusRadiodurans.fasta"},
				{"PA", root + "bacteria/GCRich/PseudomonasAeruginosaComplete.fasta"},
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
		*/
	}
}

func main() {
	sources := getSources()

	for i := 0; i < len(sources); i++ {
		source := sources[i]
		fmt.Printf("Loading %s...\n", source.name)
		g := genomes.LoadGenomes(source.fname, "", true)

		fmt.Printf("Counting %s...\n", source.name)
		for len := 1; len <= 2; len++ {
			dm := FindDinucs(g, 0, len)

			outName := fmt.Sprintf("%s-%d-nts.txt", source.name, len)
			dm.Save(outName)
			fmt.Printf("Wrote %s\n", outName)
		}
	}
}
