package genomes

import (
	"errors"
	"genomics/utils"
)

type Search interface {
	Start()
	Get() (int, error)
	IsForwards() bool	// Was the last match forwards or backwards?
	Next()
	End() bool
	GenomeLength() int
}

type LinearSearch struct {
	needle    []byte
	haystack  *Genomes
	which     int
	pos       int
	lastFound int // The last place we found a match (usually == pos-1)
	// If 0.0 then matches must be exact. Otherwise it's how many errors you
	// tolerate as proportion of needle length. So 0.1 would require a 90%
	// match.
	tolerance         float64
	allowedMismatches int

	// How many mismatches there were in the last result
	lastMismatches int
}

func (s *LinearSearch) Init(haystack *Genomes, which int, needle []byte,
	tolerance float64) {
	s.haystack = haystack
	s.which = which
	s.needle = needle
	s.allowedMismatches = int(tolerance * float64(len(needle)))
	s.Start()
}

func (s *LinearSearch) Start() {
	s.pos = 0
	s.Next()
}

func (s *LinearSearch) StartAt(pos int) {
	s.pos = pos
	s.Next()
}

func (s *LinearSearch) Get() (int, error) {
	nts := s.haystack.Nts[s.which]
	if s.lastFound < len(nts) {
		return s.lastFound, nil
	}
	return 0, errors.New("Off end")
}

func (s *LinearSearch) IsForwards() bool {
	return true
}

func (s *LinearSearch) FastForward(pos int) {
	nts := s.haystack.Nts[s.which]
	s.lastFound = len(nts)
	s.pos = pos
}

func (s *LinearSearch) Next() {
	nts := s.haystack.Nts[s.which]
	n, m := len(nts), len(s.needle)

searching:
	for ; s.pos < n-m; s.pos++ {
		var mismatches int
		for i := 0; i < m; i++ {
			if nts[s.pos+i] != s.needle[i] {
				mismatches++
				if mismatches >= s.allowedMismatches {
					continue searching
				}
			}
		}
		s.lastFound = s.pos
		s.lastMismatches = mismatches
		s.pos++
		return
	}
	s.lastFound = n
}

func (s *LinearSearch) GenomeLength() int {
	return s.haystack.Length()
}

func (s *LinearSearch) End() bool {
	nts := s.haystack.Nts[s.which]
	n, m := len(nts), len(s.needle)
	return s.pos == n-m
}

func (s *LinearSearch) Mismatches() int {
	return s.lastMismatches
}

func SearchAll(s Search) []int {
	ret := make([]int, 0)
	for ; !s.End(); s.Next() {
		pos, _ := s.Get()
		ret = append(ret, pos)
	}
	return ret
}

// Generic bidirectional search
type BidiSearch struct {
	forwards  Search
	backwards Search
	genomeLen int
}

func (s *BidiSearch) Start() {
	s.forwards.Start()
	s.backwards.Start()
}

func (s *BidiSearch) Next() {
	if !s.forwards.End() {
		s.forwards.Next()
	} else {
		s.backwards.Next()
	}
}

func (s *BidiSearch) Get() (int, error) {
	if !s.forwards.End() {
		return s.forwards.Get()
	} else {
		return s.backwards.Get()
	}
}

func (s *BidiSearch) IsForwards() bool {
	return !s.forwards.End()
}

func (s *BidiSearch) End() bool {
	return s.forwards.End() && s.backwards.End()
}

func (s *BidiSearch) GenomeLength() int {
	return s.genomeLen
}

func NewLinearSearch(haystack *Genomes,
	which int, needle []byte, tolerance float64) *LinearSearch {
	var ret LinearSearch
	ret.Init(haystack, which, needle, tolerance)
	return &ret
}

func NewIndexSearch(root string, needle []byte) *IndexSearch {
	var ret IndexSearch
	ret.Init(root, needle)
	return &ret
}

func NewBidiLinearSearch(haystack *Genomes,
	which int, needle []byte, tolerance float64) *BidiSearch {
	var forwards, backwards LinearSearch

	forwards.Init(haystack, which, needle, tolerance)
	backwards.Init(haystack, which, utils.ReverseComplement(needle), tolerance)
	return &BidiSearch{&forwards, &backwards, haystack.Length()}
}

func NewBidiIndexSearch(root string, needle []byte) *BidiSearch {
	var forwards, backwards IndexSearch

	forwards.Init(root, needle)
	backwards.Init(root, utils.ReverseComplement(needle))
	return &BidiSearch{&forwards, &backwards, forwards.GenomeLength()}
}
