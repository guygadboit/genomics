package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"io"
	"log"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

type Insertion struct {
	pos   int    // Where
	nts   []byte // What
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
		total, float64(total)/float64(len(insertions)))
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

/*
	Ponderous sort but less trouble than updating my whole system to the latest
	Go version that has slices.SortFunc
*/

type posCount struct {
	pos		int
	count	int
}

type posCounts []posCount

func (p posCounts) Len() int {
	return len(p)
}

func (p posCounts) Less(i, j int) bool {
	return p[i].pos < p[j].pos
}

func (p posCounts) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func byLocation(insertions []Insertion, minLength int) {
	positions := make(map[int]int)

	for i := 0; i < len(insertions); i++ {
		count, _ := positions[insertions[i].pos]
		positions[insertions[i].pos] = count+1
	}

	posCounts := make(posCounts, 0, len(positions))
	for k, v := range positions {
		posCounts = append(posCounts, posCount{k, v})
	}

	sort.Sort(posCounts)

	fd, err := os.Create("locations.txt")
	if err != nil {
		log.Fatal("Can't create file")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)

	for i := 0; i < posCounts.Len(); i++ {
		pc := posCounts[i]
		fmt.Fprintf(w, "%d %d\n", pc.pos, pc.count)
	}

	fmt.Printf("Wrote locations.txt\n")
}

func main() {
	insertions := LoadInsertions("insertions.txt", 9, 2)
	Summary(insertions)
	// findInVirus(insertions, 9)
	byLocation(insertions, 9)

	/*
		for i := 9; i < 50; i++ {
			findInVirus(insertions, i)
		}
	*/
}
