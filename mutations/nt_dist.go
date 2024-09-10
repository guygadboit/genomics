package mutations

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"math/rand"
)

// Counts for each nucleotide in a genome
type NucDistro struct {
	nts   map[byte]int
	total int
}

type NtIterator interface {
	Start()
	Get() byte
	Next()
	End() bool
}

type GenomeIterator struct {
	g    *genomes.Genomes
	i, j int
}

func (gi *GenomeIterator) Start() {
	gi.i, gi.j = 0, 0
}

func (gi *GenomeIterator) Get() byte {
	return gi.g.Nts[gi.i][gi.j]
}

func (gi *GenomeIterator) Next() {
	gi.j++
	if gi.j == gi.g.Length() {
		gi.j = 0
		gi.i++
	}
}

func (gi *GenomeIterator) End() bool {
	return gi.i == gi.g.NumGenomes()
}

func NewGenomeIterator(g *genomes.Genomes) NtIterator {
	ret := GenomeIterator{g, 0, 0}
	return &ret
}

func (nd *NucDistro) Count(it NtIterator) {
	for it.Start(); !it.End(); it.Next() {
		nt := it.Get()
		if !utils.IsRegularNt(nt) {
			continue
		}

		count, _ := nd.nts[nt]
		nd.nts[nt] = count + 1
		nd.total += 1
	}
}

func NewNucDistro(it NtIterator) *NucDistro {
	ret := NucDistro{nts: make(map[byte]int)}
	ret.Count(it)
	return &ret
}

func (nd *NucDistro) Show() {
	for k := range nd.nts {
		fmt.Printf("%c: %d %.2f%%\n", k, nd.nts[k],
			float64(100.0*nd.nts[k])/float64(nd.total))
	}
	fmt.Printf("Total: %d\n", nd.total)
}

/*
Pick a nucleotide randomly from the distribution represented by nd
*/
func (nd *NucDistro) Random() byte {
	r := rand.Intn(nd.total)

	var k byte
	for k = range nd.nts {
		if r < nd.nts[k] {
			break
		}
	}
	return k
}

func (nd *NucDistro) RandomSequence(s []byte) {
	for i := 0; i < len(s); i++ {
		s[i] = nd.Random()
	}
}
