package genomes

import (
	"reflect"
)

type Search struct {
	needle		[]byte
	haystack	*Genomes
	which		int
	pos			int
}

func (s *Search) Init(haystack *Genomes, which int, needle []byte) {
	s.haystack = haystack
	s.which = which
	s.needle = needle
	s.pos = 0
}

func (s *Search) Next() int {
	nts := s.haystack.nts[s.which]
	n, m := len(nts), len(s.needle)
	for ; s.pos < n - m; s.pos++ {
		if reflect.DeepEqual(nts[s.pos:s.pos+m], s.needle) {
			ret := s.pos
			s.pos++
			return ret
		}
	}
	s.pos = n
	return n
}

func (s *Search) End() bool {
	nts := s.haystack.nts[s.which]
	return s.pos >= len(nts)
}
