package genomes

import (
	"bufio"
	"fmt"
	"log"
	"os"
)

type indexData map[string][]int

type index struct {
	root  string    // Where we store the files
	data  indexData // Maps patterns to their positions
	count int       // Number of positions so far recorded
}

/*
	Write the index out to files and clear it
*/
func (index *index) save() {
	for k, v := range index.data {
		fname := fmt.Sprintf("%s/%s", index.root, k)
		fd, err := os.OpenFile(fname,
			os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0644)
		if err != nil {
			log.Fatal("Can't open file")
		}
		defer fd.Close()

		fp := bufio.NewWriter(fd)

		for i := 0; i < len(v); i++ {
			fmt.Println(fp, v[i])
		}

		fp.Flush()
		fd.Close()
	}

	index.data = make(indexData)
	index.count = 0
}

func (index *index) add(pat string, pos int) {
	_, there := index.data[pat]
	if !there {
		index.data[pat] = make([]int, 0)
	}
	index.data[pat] = append(index.data[pat], pos)
	index.count++

	// FIXME obviously this count will be bigger
	if index.count == 100 {
		index.save()
	}
}

func (index *index) Build(genome *Genomes,
	root string, prefix []byte, length int) {
	index.root = root
	index.data = make(indexData)

	nts := genome.Nts[0]
	for i := 0; i < genome.Length()-length; i++ {
		for j := 0; j < len(prefix); j++ {
			if nts[i] != prefix[j] {
				continue
			}
		}

		pat := string(nts[i : i+length])
		index.add(pat, i)
	}
}
