package genomes

import (
	"errors"
)

type Search struct {
	needle   []byte
	haystack *Genomes
	which    int
	pos      int
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
	if s.pos < len(nts) {
		return s.pos, nil
	}
	return 0, errors.New("Off end")
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
		s.pos++
		return
	}
}

func (s *Search) End() bool {
	nts := s.haystack.Nts[s.which]
	n, m := len(nts), len(s.needle)
	return s.pos == n-m
}
