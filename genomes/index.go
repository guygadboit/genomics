package genomes

import (
	"bufio"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"strings"
	"errors"
)

type indexData map[string][]int

type Index struct {
	root    string    // Where we store the files
	data    indexData // Maps patterns to their positions
	wordLen int       // The word length we're using (usually 6)
	count   int       // Number of positions so far recorded
}

/*
	Write the index out to files and clear it
*/
func (index *Index) Save() {
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
		index.Save()
	}
}

func (index *Index) Build(genome *Genomes,
	root string, length int, verbose bool) {
	index.root = root
	index.data = make(indexData)
	index.wordLen = length
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

func (index *Index) Init(root string, wordLen int) {
	index.root = root
	index.wordLen = wordLen
}

type IndexSearch struct {
	haystack *Genomes
	which    int
	needle   []byte
	index    *Index

	// The places where the first wordLen nts are found, then the second
	// wordLen nts, etc.
	positions [][]int

	// Where we are in our iteration through positions
	iterators []int

	// Search the last few nts (that are less than a whole wordLen)
	atLastMile bool
	lastMile   Search
}

/*
	Read in one of our cache files. This contains all the positions of a
	particular pattern, like CCGGGT or whatever. Return the positions.
*/
func readFile(root string, pattern string) []int {
	ret := make([]int, 0)
	f := utils.NewFileReader(fmt.Sprintf("%s/%s", root, pattern))

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

func (s *IndexSearch) Init(root string, wordLen int,
	haystack *Genomes, which int, needle []byte) {
	s.haystack = haystack
	s.which = which
	s.needle = needle
	s.index = &Index{root: root, wordLen: wordLen}

	m := s.index.wordLen

	// Depth is how deep we can go with the cached positions-- if wordLen is 6
	// say and we're searching for something 13 long, we will go into the cache
	// for the first 6 nts, then again for the next 6, and then search the last
	// 1 with our "last mile" search. So in that case depth will be 2.
	depth := len(needle) / m
	s.positions = make([][]int, depth)

	for i := 0; i < depth; i++ {
		s.positions[i] = readFile(s.index.root, string(needle[i*m:i*m+m]))
	}

	s.Start()
}

func (s *IndexSearch) Start() {
	s.iterators = make([]int, len(s.positions))
	s.Next()
}

/*
	FIXME: OK you can make this do binary searches to increment to the next
	position that might work. The idea is that whenever you increase a counter
	it needs to have a chance of being +6 from the previous one.
*/
func (s *IndexSearch) incr() bool {
	for i := 0; i < len(s.iterators); i++ {
		s.iterators[i]++
		if s.iterators[i] < len(s.positions[i]) {
			return true
		}
		for j := 0; j <= i; j++ {
			s.iterators[j] = 0
		}
	}
	return false
}

/*
	Do the current iterators record a match? If so record the position. If not
	then -1.
*/
func (s *IndexSearch) haveMatch() int {
	m := s.index.wordLen

	// This is the case where the needle is shorter than wordLen
	if len(s.iterators) == 0 {
		return 0
	}

	for i := 0; i < len(s.iterators); i++ {
		if i > 1 {
			if s.positions[i][s.iterators[i]] !=
				s.positions[i-1][s.iterators[i-1]]+m {
				return -1
			}
		}
	}

	return s.positions[0][s.iterators[0]]
}

func (s *IndexSearch) Next() {
	for {
		if s.atLastMile {
			s.lastMile.Next()
			if s.lastMile.End() {
				s.atLastMile = false
			}
			return
		}

		pos := s.haveMatch()
		if pos != -1 {
			s.lastMile.Init(s.haystack, s.which, s.needle, 0.0)
			s.lastMile.FastForward(pos)
			s.atLastMile = true
			s.Next()
			return
		}
		if !s.incr() {
			break
		}
	}
}

func (s *IndexSearch) Get() (int, error) {
	if s.atLastMile {
		return s.lastMile.Get()
	}
	return 0, errors.New("Off end")
}

func (s *IndexSearch) End() bool {
	if s.atLastMile {
		return s.lastMile.End()
	}
	return true
}
