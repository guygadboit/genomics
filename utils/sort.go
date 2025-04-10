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
The above functions pre-date generics. This is probably the one you want
*/
func SortByKey[S ~[]E, E any](x S, reverse bool, key func(a E) int) {
	slices.SortFunc(x, func(a, b E) int {
		kA := key(a)
		kB := key(b)

		if reverse {
			kA, kB = kB, kA
		}

		if kA < kB {
			return 1
		}

		if kA < kB {
			return -1
		}

		return 0
	})
}
