package degeneracy

import (
	"genomics/genomes"
)

/*
For each codon, the "foldness of its degeneracy": in other words, how many
alternative nts in the 3rd position would give you the same AA?
*/
var DegeneracyTable map[string]int

type Codon struct {
	genomes.Codon
	Fold int
}

type Translation []Codon

func newDegeneracyTable() map[string]int {
	ret := make(map[string]int)
	for k, v := range genomes.CodonTable {
		for _, syn := range genomes.ReverseCodonTable[v] {
			if syn[:2] == k[:2] {
				ret[k]++
			}
		}
	}
	return ret
}

func init() {
	DegeneracyTable = newDegeneracyTable()
}

func AddDegeneracy(t genomes.Translation) Translation {
	ret := make([]Codon, len(t))
	for i, c := range t {
		ret[i] = Codon{c, DegeneracyTable[c.Nts]}
	}
	return ret
}
