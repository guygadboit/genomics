package genomes

import (
	"errors"
)

type SearchIf interface {
	Start()
	Get() (int, error)
	Next()
	End() bool
}

type Search struct {
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
}

func (s *Search) Init(haystack *Genomes, which int, needle []byte,
	tolerance float64) {
	s.haystack = haystack
	s.which = which
	s.needle = needle
	s.allowedMismatches = int(tolerance * float64(len(needle)))
	s.Start()
}

func (s *Search) Start() {
	s.pos = 0
	s.Next()
}

func (s *Search) Get() (int, error) {
	nts := s.haystack.Nts[s.which]
	if s.lastFound < len(nts) {
		return s.lastFound, nil
	}
	return 0, errors.New("Off end")
}

func (s *Search) FastForward(pos int) {
	nts := s.haystack.Nts[s.which]
	s.lastFound = len(nts)
	s.pos = pos
}

func (s *Search) Next() {
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
		s.pos++
		return
	}
	s.lastFound = n
}

func (s *Search) End() bool {
	nts := s.haystack.Nts[s.which]
	n, m := len(nts), len(s.needle)
	return s.pos == n-m
}

func SearchAll(s SearchIf) []int {
	ret := make([]int, 0)
	for ; !s.End(); s.Next() {
		pos, _ := s.Get()
		ret = append(ret, pos)
	}
	return ret
}
