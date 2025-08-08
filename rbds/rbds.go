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

type NtComparison struct {
	Comparison
	N, S int
}

func (c *NtComparison) SNRatio() float64 {
	return float64(c.S) / float64(c.N)
}

func PrintComparisons(c []Comparison, g *genomes.Genomes, caption string) {
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

func PrintNtComparisons(c []NtComparison,
	g *genomes.Genomes, caption string) NtComparison {
	slices.SortFunc(c, func(a, b NtComparison) int {
		ar, br := a.SNRatio(), b.SNRatio()
		if ar < br {
			return 1
		}
		if ar > br {
			return -1
		}
		return 0
	})
	for _, c := range c {
		fmt.Printf("%s %s %.4f S=%d N=%d %.4f\n", caption,
			g.Names[c.index], c.similarity, c.S, c.N, c.SNRatio())
	}
	return c[0]
}

/*
orfStart is the start of the spike; start and are the offsets into the spike
of the bit you're interested in (which will be the RBD or the RBM)
*/
func Compare(g *genomes.Genomes, ref int, orfStart utils.OneBasedPos,
	start, end utils.OneBasedPos, protein bool) {
	comparisons := make([]Comparison, 0)
	s := int((orfStart + start*3) - 1)
	e := int(orfStart + end*3)
	fmt.Println(s, e)
	g.Truncate(s, e)
	g.ResetOrfs()

	fmt.Printf("Ref: %d\n", ref)

	// Protein first
	if protein {
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
		PrintComparisons(comparisons, g, "Protein")
	}

	// g.SaveMulti("RBDs.fasta")
	ntComparisons := make([]NtComparison, 0)
	for i := 0; i < g.NumGenomes(); i++ {
		if i == ref {
			continue
		}
		same := 0
		S, N := 0, 0
		for j := 0; j < g.Length(); j++ {
			if g.Nts[ref][j] == g.Nts[i][j] {
				same++
			} else {
				if silent, _, _ := genomes.IsSilent(g, j, 1, ref, i); silent {
					S++
				} else {
					N++
				}
			}
		}
		similarity := float64(same) / float64(g.Length())
		ntComparisons = append(ntComparisons,
			NtComparison{Comparison{i, similarity}, N, S})
	}
	highest := PrintNtComparisons(ntComparisons, g, "Nucleotide")
	fmt.Printf("%.4f\t%d\t%d\t%s vs %s G\n",
		highest.SNRatio(), highest.S, highest.N, g.Names[ref], g.Names[highest.index])
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

func MakeChimera(g *genomes.Genomes) *genomes.Genomes {
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

	for ref := 0; ref < g.NumGenomes(); ref++ {
		/*
		fmt.Println("RBD")
		Compare(g.Clone(), ref, SPIKE_START, rbd[0], rbd[1], false)
		*/

		fmt.Println("RBM")
		Compare(g.Clone(), ref, SPIKE_START, rbm[0], rbm[1], false)
	}

	return

	/*
		fmt.Println("ORF8")
		Compare(g.Clone(), ORF8_START, 1, 122, false)

		fmt.Println("M")
		Compare(g.Clone(), M_START, 1, 223, false)

		fmt.Println("ORF10")
		Compare(g.Clone(), ORF10_START, 1, 39, false)
		fmt.Println("N")
		Compare(g.Clone(), N_START, 1, 420, false)
	*/

	/*
		These guys are all quite a bit different from SARS2 and RaTG13 it looks
		like. But simlar to each other. They're generally similar to SARS1 than
		SARS2. Particularly Rs5725
	*/
	g = genomes.LoadGenomes("/fs/f/genomes/viruses/"+
		"suppressed_genomes/S/WH1-alignment.fasta",
		"", false)

	fmt.Println("RBD (additional)")
	Compare(g.Clone(), 0, 0, rbd[0], rbd[1], false)
	fmt.Println("RBM (additional)")
	Compare(g.Clone(), 0, 0, rbm[0], rbm[1], false)
}
