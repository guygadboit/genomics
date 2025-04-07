package utils

import (
	"sort"
	"slices"
)

type sorter struct {
	len		int
	reverse bool
	swap	func(int, int)

	// Must supply one or the other of these
	lt      func(int, int) bool
	key     func(int) float64
}

func (s *sorter) init(len int, reverse bool, swap func(int, int),
	lt func(int, int) bool, key func(int) float64) {

	s.len = len
	s.reverse = reverse
	s.swap = swap

	if lt != nil {
		s.lt = lt
	} else {
		s.key = key
	}
}

func (s *sorter) Len() int {
	return s.len
}

func (s *sorter) Less(i, j int) bool {
	var result bool

	if s.lt != nil {
		result = s.lt(i, j)
	} else {
		a, b := s.key(i), s.key(j)
		result = a < b
	}

	if s.reverse {
		return !result
	} else {
		return result
	}
}

func (s *sorter) Swap(i, j int) {
	s.swap(i, j)
}

/*
	Sort a slice. Either of lt or key can be nil. This pre-dated generics is
	probably useless now.
*/
func LegacySort(len int, reverse bool, swap func(int, int),
	lt func(int, int) bool, key func(int) float64) {
	var s sorter
	s.init(len, reverse, swap, lt, key)
	sort.Sort(&s)
}

/*
The above pre-dated generics. This is the Sort function you actually probably
want
*/
func Sort[T any](s []T, reverse bool, key func(t T) int) {
	var ret int
	if reverse {
		ret = -1
	} else {
		ret = 1
	}
	slices.SortFunc(s, func(a, b T) int {
		if key(a) > key(b) {
			return ret
		}
		if key(a) < key(b) {
			return -ret
		}
		return 0
	})
}
