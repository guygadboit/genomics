package main

import (
	"os"
	"log"
	"bufio"
	"strconv"
	"fmt"
	"io"
	"strings"
	"regexp"
	"genomics/genomes"
)

type Insertion struct {
	pos  int    // Where
	nts  []byte // What
	nSeqs int    // How many times
}

func (i *Insertion) ToString() string {
	return fmt.Sprintf("%s at %d (%d seqs)", string(i.nts),
		i.pos, i.nSeqs)
}

func atoi(s string) int {
	ret, err := strconv.Atoi(s)
	if err != nil {
		log.Fatal("Bad integer")
	}
	return ret
}

func LoadInsertions(fname string, minLen int, minSeqs int) []Insertion {
	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal("Can't open %s", fname)
	}
	defer fd.Close()

	ret := make([]Insertion, 0)
	fp := bufio.NewReader(fd)

	pat := regexp.MustCompile(`ins_(\d+):([A-Z]+) \((\d+) seqs\)`)

reading:
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break reading
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)
		groups := pat.FindAllStringSubmatch(line, -1)

		ins := Insertion{atoi(groups[0][1]),
			[]byte(groups[0][2]),
			atoi(groups[0][3])}

		if len(ins.nts) < minLen {
			continue
		}

		if ins.nSeqs < minSeqs {
			continue
		}

		ret = append(ret, ins)
	}

	return ret
}

func Summary(insertions []Insertion) {
	fmt.Printf("%d insertions\n", len(insertions))

	var total int
	for i := 0; i < len(insertions); i++ {
		total += len(insertions[i].nts)
	}
	fmt.Printf("%d nts altogether (average length %.2f)\n",
		total, float64(total) / float64(len(insertions)))
}

func findInVirus(insertions []Insertion, minLength int) {
	wh1 := genomes.LoadGenomes("../fasta/WH1.fasta",
		"../fasta/WH1.orfs")
	var search genomes.Search

	var found, count int
searching:
	for i := 0; i < len(insertions); i++ {
		nts := insertions[i].nts
		if len(nts) < minLength {
			continue
		}
		count++

		for search.Init(wh1, 0, nts); !search.End(); search.Next() {
			found++
			continue searching
		}

		rc := genomes.ReverseComplement(insertions[i].nts)
		for search.Init(wh1, 0, rc); !search.End(); search.Next() {
			found++
			continue searching
		}

		ins := insertions[i]
		fmt.Printf("ins_%d:%s (%d seqs)\n", ins.pos, ins.nts, ins.nSeqs)
	}

	fmt.Printf("Length %d: %d (/%d) were found in SC2 itself\n", minLength,
		found, count)
}

func main() {
	insertions := LoadInsertions("insertions.txt", 9, 2)
	Summary(insertions)
	findInVirus(insertions, 9)

	/*
	for i := 9; i < 50; i++ {
		findInVirus(insertions, i)
	}
	*/
}
