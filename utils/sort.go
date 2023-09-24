package utils

import (
	"sort"
)

type sorter struct {
	data []interface{}
	lt   func(int, int) bool
	key  func(int) float64
}

func (s sorter) init(data []interface{},
	lt func(int, int) bool,
	key func(int) float64) {

	s.data = data
	if lt != nil {
		s.lt = lt
	} else {
		s.key = key
	}
}

func (s sorter) Len() int {
	return len(s.data)
}

func (s sorter) Less(i, j int) bool {
	if s.lt != nil {
		return s.lt(i, j)
	} else {
		a, b := s.key(i), s.key(j)
		return a < b
	}
}

func (s sorter) Swap(i, j int) {
	s.data[i], s.data[j] = s.data[j], s.data[i]
}

/*
	Sort a slice. Either of lt or key can be nil
*/
func Sort(data []interface{}, lt func(int, int) bool, key func(int) float64) {
	var s sorter
	s.init(data, lt, key)
	sort.Sort(s)
}
