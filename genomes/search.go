package genomes

import (
	"errors"
	"reflect"
)

type Search struct {
	needle   []byte
	haystack *Genomes
	which    int
	pos      int
}

func (s *Search) Init(haystack *Genomes, which int, needle []byte) {
	s.haystack = haystack
	s.which = which
	s.needle = needle
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
	for ; s.pos < n-m; s.pos++ {
		if reflect.DeepEqual(nts[s.pos:s.pos+m], s.needle) {
			return
		}
	}
}

func (s *Search) End() bool {
	nts := s.haystack.Nts[s.which]
	n, m := len(nts), len(s.needle)
	return s.pos == n-m
}
