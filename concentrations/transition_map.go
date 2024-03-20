package main

import (
	"bufio"
	"fmt"
	"genomics/stats"
	"slices"
)

// A->B
type DirectedTransition struct {
	ANts string
	BNts string
}

/*
Convert a directed transition to an undirected one, rearranging its ANts and
BNts into a "canonical order" so that either way round will give you the same
map key.
*/
func (t DirectedTransition) Undirected() UndirectedTransition {
	var a, b string
	if t.ANts < t.BNts {
		a, b = t.ANts, t.BNts
	} else {
		a, b = t.BNts, t.ANts
	}
	return UndirectedTransition{a, b}
}

// A<->B which means either A->B or B->A, we don't know which
type UndirectedTransition DirectedTransition

func (t DirectedTransition) String() string {
	return fmt.Sprintf("%s->%s", t.ANts, t.BNts)
}

func (t UndirectedTransition) String() string {
	return fmt.Sprintf("%s<->%s", t.ANts, t.BNts)
}

type Transition interface {
	DirectedTransition | UndirectedTransition
	String() string
}

type TransitionCounter[T Transition] struct {
	Counts map[T]int
	Total  int

	// The keys into Counts sorted by most frequent first
	OrderedKeys []T
}

func (t *TransitionCounter[T]) Init() {
	t.Counts = make(map[T]int)
}

func (t *TransitionCounter[T]) updateTotal() {
	t.Total = 0
	for _, v := range t.Counts {
		t.Total += v
	}
}

func (t *TransitionCounter[T]) orderKeys() {
	keys := make([]T, 0)
	for k, _ := range t.Counts {
		keys = append(keys, k)
	}
	slices.SortFunc(keys, func(a, b T) int {
		return t.Counts[b] - t.Counts[a]
	})
	t.OrderedKeys = keys
}

func (t *TransitionCounter[T]) finalize() {
	t.updateTotal()
	t.orderKeys()
}

type TransitionMap struct {
	Directed   TransitionCounter[DirectedTransition]
	Undirected TransitionCounter[UndirectedTransition]
}

func (tm *TransitionMap) Init() {
	tm.Directed.Init()
	tm.Undirected.Init()
}

func (tm *TransitionMap) Add(aNts, bNts string) {
	d := DirectedTransition{aNts, bNts}
	u := d.Undirected()
	tm.Directed.Counts[d]++
	tm.Undirected.Counts[u]++
}

func (tm *TransitionMap) Combine(other *TransitionMap) {
	for k, v := range other.Directed.Counts {
		tm.Directed.Counts[k] += v
	}
	for k, v := range other.Undirected.Counts {
		tm.Undirected.Counts[k] += v
	}
}

// Call this when you've finished adding things to the map
func (tm *TransitionMap) Finalize() {
	tm.Directed.finalize()
	tm.Undirected.finalize()
}

/*
Given two transition maps (a real one and simulated one typically) return the
significance of how the real one differs from the simulated one.
*/
func CompareTransition[T Transition](
	a, b *TransitionCounter[T], t T) stats.ContingencyTable {
	A, B := a.Counts[t], a.Total-a.Counts[t]
	C, D := b.Counts[t], b.Total-b.Counts[t]

	var ret stats.ContingencyTable
	ret.Init(A, B, C, D)
	return ret
}

func compareTransitionCounters[T Transition](
	a, b *TransitionCounter[T], w *bufio.Writer) {
	for i, k := range a.OrderedKeys {
		fmt.Fprintf(w, "%s %d", k, a.Counts[k])
		if i < 6 {
			ct := CompareTransition(a, b, k)
			ct.FisherExact()
			fmt.Fprintf(w, " OR=%.4f p=%g\n", ct.OR, ct.P)
		} else {
			fmt.Fprintf(w, "\n")
		}
	}
}

// Report on a, and compare it to b
func CompareTransitionMaps(a, b *TransitionMap, w *bufio.Writer) {
	fmt.Fprintln(w, "Undirected transitions")
	compareTransitionCounters(&a.Undirected, &b.Undirected, w)

	fmt.Fprintln(w, "\nDirected transitions")
	compareTransitionCounters(&a.Directed, &b.Directed, w)
}

type GraphDatum struct {
	Key  string
	Real float64
	Sim  float64
}

func GraphData(realMap, simMap *TransitionMap, w *bufio.Writer) {
	data := make([]GraphDatum, 0)

	for k, v := range realMap.Undirected.Counts {
		r := float64(v) / float64(realMap.Undirected.Total)

		simV := simMap.Undirected.Counts[k]
		s := float64(simV) / float64(simMap.Undirected.Total)

		data = append(data, GraphDatum{k.String(), r, s})
	}

	// Sort these alphabetically, so we can compare graphs from different
	// viruses more easily.
	slices.SortFunc(data, func(a, b GraphDatum) int {
		if a.Key < b.Key {
			return -1
		}
		return 1
	})

	for _, d := range data {
		fmt.Fprintf(w, "%s %.4f %.4f\n", d.Key, d.Real, d.Sim)
	}
}
