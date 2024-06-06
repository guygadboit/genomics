package main

import (
	"fmt"
	"genomics/database"
)

type Transition struct {
	From byte
	To   byte
}

type TransitionCounter struct {
	counter map[Transition]int
	seen    map[database.Mutation]bool
}

func NewTransitionCounter() *TransitionCounter {
	var ret TransitionCounter
	ret.counter = make(map[Transition]int)
	ret.seen = make(map[database.Mutation]bool)
	return &ret

}

func (tc *TransitionCounter) Add(mut database.Mutation) {
	seen := tc.seen[mut]
	if seen {
		return
	}

	t := Transition{mut.From, mut.To}
	tc.counter[t] += 1
	tc.seen[mut] = true
}

func (tc TransitionCounter) Print() {
	var total float64
	for _, count := range tc.counter {
		total += float64(count)
	}

	for t, count := range tc.counter {
		freq := float64(count) / total
		fmt.Printf("%c->%c: %.2f\n", t.From, t.To, freq)
	}
}
