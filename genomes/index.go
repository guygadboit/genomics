package genomes

import (
	"bufio"
	"errors"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"sort"
	"strings"
)

type indexData map[string][]int

type Index struct {
	root      string    // Where we store the files
	data      indexData // Maps patterns to their positions
	wordLen   int       // The word length we're using (usually 6)
	count     int       // Number of positions so far recorded
	genomeLen int
}

func (index *Index) saveMetadata() {
	fname := fmt.Sprintf("%s/info", index.root)
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal("Can't write index metadata")
	}
	defer fd.Close()

	w := bufio.NewWriter(fd)
	fmt.Fprintf(w, "genomeLen: %d wordLen: %d\n",
		index.genomeLen, index.wordLen)
	w.Flush()
}

func (index *Index) loadMetadata() {
	fname := fmt.Sprintf("%s/info", index.root)
	f := utils.NewFileReader(fname)
	defer f.Close()

	line, err := f.ReadString('\n')
	if err != nil {
		log.Fatal("Can't read index metadata")
	}

	line = strings.TrimSpace(line)
	fields := strings.Fields(line)
	index.genomeLen = utils.Atoi(fields[1])
	index.wordLen = utils.Atoi(fields[3])
}

/*
	Write the index out to files and clear it
*/
func (index *Index) save() {
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
			fmt.Fprintln(fp, v[i])
		}

		fp.Flush()
		fd.Close()
	}

	index.data = make(indexData)
	index.count = 0
}

func (index *Index) add(pat string, pos int) {
	_, there := index.data[pat]
	if !there {
		index.data[pat] = make([]int, 0)
	}
	index.data[pat] = append(index.data[pat], pos)
	index.count++

	if index.count == 1024*1024 {
		index.save()
	}
}

func (index *Index) Build(genome *Genomes,
	root string, length int, verbose bool) {
	index.root = root
	index.data = make(indexData)
	index.wordLen = length
	index.genomeLen = genome.Length()
	n := genome.Length() - length

	nts := genome.Nts[0]
	for i := 0; i < n; i++ {
		pat := string(nts[i : i+length])
		index.add(pat, i)

		if verbose && i%10000 == 0 {
			fmt.Printf("%.2f%%\n", float64(i*100)/float64(n))
		}
	}
}

func (index *Index) Save() {
	index.save()
	index.saveMetadata()
}

func (index *Index) Init(root string, wordLen int) {
	index.root = root
	index.wordLen = wordLen
}

type IndexSearch struct {
	needle []byte
	index  *Index

	// The places where the first wordLen nts are found, then the second
	// wordLen nts, etc.
	positions [][]int

	// Where we are in our iteration through the first list of positions
	iterator int

	// The offset to use for the last set of positions (the others are all
	// +wordLen successively from the start). This offset is relative to
	// +wordLen so that the default value of 0 does nothing.
	lastOffset int

	// Where we last found something (usually == pos-1)
	lastFound int
}

/*
	Read in one of our cache files. This contains all the positions of a
	particular pattern, like CCGGGT or whatever. Return the positions.
*/
func readFile(root string, pattern string) []int {
	ret := make([]int, 0)

	fname := fmt.Sprintf("%s/%s", root, pattern)
	if _, err := os.Stat(fname); err != nil {
		fname += ".gz"
	}

	f := utils.NewFileReader(fname)
	defer f.Close()

loop:
	for {
		line, err := f.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)
		ret = append(ret, utils.Atoi(line))
	}
	return ret
}

func (s *IndexSearch) Init(root string, needle []byte) {
	s.needle = needle
	n := len(needle)

	s.index = &Index{root: root}
	s.index.loadMetadata()

	m := s.index.wordLen

	depth := n / m

	// Capacity of +1 because if there is an overhang we will add another list
	// of positions below.
	s.positions = make([][]int, depth, depth+1)

	for i := 0; i < depth; i++ {
		s.positions[i] = readFile(s.index.root, string(needle[i*m:i*m+m]))
	}

	/*
		If the needle is not a multiple of the index wordLen we need an
		additional fragment to match, which is the *last* wordLen nts of the
		needle. This overlaps with the previous thing we looked for. But that's
		fine. We just have to correct the target offset we will be looking for
		for these bits. We save that correction in lastOffset.
	*/
	overhang := n - depth*m
	if overhang > 0 {
		s.lastOffset = -m + overhang
		s.positions = append(s.positions,
			readFile(s.index.root, string(needle[n-m:n])))
	}

	s.Start()
}

func (s *IndexSearch) Start() {
	s.lastFound = -1
	s.iterator = 0
	s.Next()
}

func (s *IndexSearch) incr() bool {
	nWords := len(s.positions)

	for ; s.iterator < len(s.positions[0]); {
		prev := s.positions[0][s.iterator]
		found := true

		for i := 1; i < nWords; i++ {
			target := prev + s.index.wordLen
			if i == nWords-1 {
				target += s.lastOffset
			}

			n := len(s.positions[i])
			pos := sort.Search(n, func(j int) bool {
				return s.positions[i][j] >= target
			})
			if pos == n {
				found = false
				break
			}
			if s.positions[i][pos] != target {
				found = false
				break
			}
			prev = s.positions[i][pos]
		}

		if found {
			s.lastFound = s.iterator
			s.iterator++
			return true
		} else {
			s.iterator++
		}
	}
	return false
}

func (s *IndexSearch) Next() {
	s.incr()
}

func (s *IndexSearch) Get() (int, error) {
	if s.lastFound != -1 {
		return s.positions[0][s.lastFound], nil
	}
	return 0, errors.New("Off end")
}

func (s *IndexSearch) End() bool {
	return s.iterator == len(s.positions[0])
}

func (s *IndexSearch) GenomeLength() int {
	return s.index.genomeLen
}
