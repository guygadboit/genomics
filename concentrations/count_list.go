package main

import (
	"fmt"
	"slices"
)

type Count[T comparable] struct {
	Key   T
	Count int
}

type CountList[T comparable] []Count[T]

/*
Turn a map of anything to int counts into a CountList, just so that we can
easily sort it and then print it out.
*/
func NewCountList[T comparable](m map[T]int) CountList[T] {
	ret := make(CountList[T], 0)
	for k, v := range m {
		ret = append(ret, Count[T]{k, v})
	}
	return ret
}

func SortCountList[T comparable](cl CountList[T]) {
	slices.SortFunc(cl, func(a, b Count[T]) int {
		return b.Count - a.Count
	})
}

func SortPrintCountList[T comparable](cl CountList[T]) {
	SortCountList(cl)
	for _, c := range cl {
		fmt.Printf("%s  %d\n", c.Key, c.Count)
	}
}
