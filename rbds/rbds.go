package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"slices"
)

var SPIKE_START utils.OneBasedPos = utils.OneBasedPos(21563)
var ORF8_START utils.OneBasedPos = utils.OneBasedPos(27894)
var M_START utils.OneBasedPos = utils.OneBasedPos(26523)
var ORF10_START utils.OneBasedPos = utils.OneBasedPos(29558)
var N_START utils.OneBasedPos = utils.OneBasedPos(28274)

type Comparison struct {
	index      int
	similarity float64
}

func printComparisons(g *genomes.Genomes, c []Comparison, caption string) {
	slices.SortFunc(c, func(a, b Comparison) int {
		if a.similarity < b.similarity {
			return 1
		}
		if a.similarity > b.similarity {
			return -1
		}
		return 0
	})
	for _, c := range c {
		fmt.Printf("%s %s %.4f\n", caption, g.Names[c.index], c.similarity)
	}
}

/*
orfStart is the start of the spike; start and are the offsets into the spike
of the bit you're interested in (which will be the RBD or the RBM)
*/
func Compare(g *genomes.Genomes, orfStart utils.OneBasedPos,
	start, end utils.OneBasedPos) {
	comparisons := make([]Comparison, 0)
	s := int((orfStart + start*3) - 1)
	e := int(orfStart + end*3)
	fmt.Println(s, e)
	g.Truncate(s, e)
	g.ResetOrfs()

	// Protein first
	ref := 0
	refTrans := genomes.Translate(g, ref)

	for i := 0; i < g.NumGenomes(); i++ {
		if i == ref {
			continue
		}
		same := 0
		trans := genomes.Translate(g, i)
		for j, codon := range refTrans {
			if codon.Aa == trans[j].Aa {
				same++
			}
		}
		similarity := float64(same) / float64(len(refTrans))
		comparisons = append(comparisons, Comparison{i, similarity})
	}
	printComparisons(g, comparisons, "Protein")

    // g.SaveMulti("RBDs.fasta")
	comparisons = make([]Comparison, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		if i == ref {
			continue
		}
		same := 0
		for j := 0; j < g.Length(); j++ {
			if g.Nts[ref][j] == g.Nts[i][j] {
				same++
			}
		}
		similarity := float64(same) / float64(g.Length())
		comparisons = append(comparisons, Comparison{i, similarity})
	}
	printComparisons(g, comparisons, "Nucleotide")
}

func swap() {
	g := genomes.LoadGenomes("rbm.fasta", "", false)
	g = g.Swap(0, 8)
	g.SaveMulti("rbm-ratfirst.fasta")

	g = genomes.LoadGenomes("rbd.fasta", "", false)
	g = g.Swap(0, 8)
	g.SaveMulti("rbd-ratfirst.fasta")
}

func centroid() {
	g := genomes.LoadGenomes("rbm.fasta", "", false)
	c := g.AACentroid()
	for i := 0; i < g.NumGenomes(); i++ {
		d := utils.VecDistance(g.ToAAVector(i), c)
		fmt.Printf("%d: %.4f\n", i, d)
	}
}

func MakeChimera(g* genomes.Genomes) *genomes.Genomes {
	ret := g.Filter(0, 7, 8, 33)
	ret.DeepCopy(2)
	for i := 22492; i < 23117; i++ {
		ret.Nts[2][i] = ret.Nts[3][i]
	}
	ret.Nts = ret.Nts[0:3]
	ret.Names = ret.Names[0:3]
	ret.Names[2] = "Chimera"
	ret.SaveMulti("chimera.fasta")
	return ret
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives-short-names.fasta",
		"../fasta/WH1.orfs", false)

	MakeChimera(g)

	// As AA offsets into S
	rbd := []utils.OneBasedPos{333, 679}
	rbm := []utils.OneBasedPos{438, 506}

	fmt.Println("RBD")
	Compare(g.Clone(), SPIKE_START, rbd[0], rbd[1])
	fmt.Println("RBM")
	Compare(g.Clone(), SPIKE_START, rbm[0], rbm[1])

	fmt.Println("ORF8")
	Compare(g.Clone(), ORF8_START, 1, 122)

	fmt.Println("M")
	Compare(g.Clone(), M_START, 1, 223)

	fmt.Println("ORF10")
	Compare(g.Clone(), ORF10_START, 1, 39)
	fmt.Println("N")
	Compare(g.Clone(), N_START, 1, 420)
    return

	/*
	These guys are all quite a bit different from SARS2 and RaTG13 it looks
	like. But simlar to each other. They're generally similar to SARS1 than
	SARS2. Particularly Rs5725
	*/
	g = genomes.LoadGenomes("/fs/f/genomes/viruses/"+
		"suppressed_genomes/spikes/WH1-alignment.fasta",
		"", false)

	fmt.Println("RBD (additional)")
	Compare(g.Clone(), 0, rbd[0], rbd[1])
	fmt.Println("RBM (additional)")
	Compare(g.Clone(), 0, rbm[0], rbm[1])
}
