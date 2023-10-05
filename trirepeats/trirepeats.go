package main

import (
	"fmt"
	"genomics/genomes"
	"log"
	"os"
	"bufio"
)

type TriRepeatMap map[string]int

func (m TriRepeatMap) Init() {
	nts := []byte{'G', 'A', 'T', 'C'}
	pat := make([]byte, 3)

	for i := 0; i < len(nts); i++ {
		for j := 0; j < len(nts); j++ {
			for k := 0; k < len(nts); k++ {
				pat[0] = nts[i]
				pat[1] = nts[j]
				pat[2] = nts[k]
				m[string(pat)] = 0
			}
		}
	}
}

func (m TriRepeatMap) Save(fname string) {
	f, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer f.Close()

	w := bufio.NewWriter(f)

	for k, v := range(m) {
		fmt.Fprintf(w, "%s: %d\n", k, v)
	}

	w.Flush()
}

func CountTriRepeats(genome *genomes.Genomes) TriRepeatMap {
	ret := make(TriRepeatMap)
	ret.Init()

	for i := 0; i < genome.Length() - 6; i++ {
		pat := string(genome.Nts[0][i:i+6])
		if pat[:3] != pat[3:] {
			continue
		}
		triNt := pat[:3]
		count := ret[triNt]
		ret[triNt] = count+1
	}
	return ret
}

func main() {
	g := genomes.LoadGenomes(os.Args[1], "", true)
	tm := CountTriRepeats(g)
	tm.Save(os.Args[2])
}

