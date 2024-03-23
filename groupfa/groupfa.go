package main

import (
	"fmt"
	"log"
	"flag"
	"slices"
	"sort"
	"strings"
	"genomics/genomes"
)

type SimilarityMatrix [][]float64

func (m SimilarityMatrix) Get(i, j int) float64 {
	// It's a symmetric matrix with 1.0 on the diagonal
	if i == j {
		return 1.0
	}
	if j < i {
		return m[i][j]
	}
	return m[j][i]
}

func (m SimilarityMatrix) Dim() int {
	return len(m)
}

func (m SimilarityMatrix) Print() {
	for i := 0; i < m.Dim(); i++ {
		for j := 0; j < m.Dim(); j++ {
			fmt.Printf("%.2f ", m.Get(i, j))
		}
		fmt.Printf("\n")
	}
}

// Make a matrix of comparisons
func CompareGenomes(g *genomes.Genomes) SimilarityMatrix {
	ret := make(SimilarityMatrix, g.NumGenomes())
	for i := 1; i < g.NumGenomes(); i++ {
		ret[i] = make([]float64, i)
		for j := 0; j < i; j++ {
			ret[i][j] = g.SequenceSimilarity(i, j)
		}
	}
	return ret
}

type Row struct {
	index	int
	total	float64
}

// Find numGroups indices that are as far apart as possible
func findBasis(m SimilarityMatrix, numGroups int) []int {
	if numGroups >= m.Dim() {
		log.Fatal("Too many groups")
	}

	rows := make([]Row, m.Dim())
	ret := make([]int, numGroups)

	for i := 0; i < m.Dim(); i++ {
		var total float64
		for j := 0; j < m.Dim(); j++ {
			total += m.Get(i, j)
		}
		rows[i] = Row{i, total}
	}

	slices.SortFunc(rows, func(a, b Row) int {
		if a.total < b.total {
			return -1
		}
		return 1
	})

	fmt.Println(rows)
	for i := 0; i < len(ret); i++ {
		ret[i] = rows[i].index
	}
	return ret
}

type FriendSet map[int]bool

func findFriends(m SimilarityMatrix, rowNum int,
	friends FriendSet, threshold float64) {
	if friends[rowNum] {
		return
	}
	friends[rowNum] = true
	for i := 0; i < m.Dim(); i++ {
		if i == rowNum {
			continue
		}
		if m.Get(rowNum, i) >= threshold {
			friends[i] = true
			findFriends(m, i, friends, threshold)
		}
	}
}

func (f FriendSet) Union(other FriendSet) {
	for k, _ := range other {
		f[k] = true
	}
}

/*
Cluster genomes into groups defined by the members of each group being
close to each other
*/
func GroupGenomes(g *genomes.Genomes, threshold float64) [][]int {
	m := CompareGenomes(g)
	visited := make(FriendSet)
	groups := make([][]int, 0)

	for i := 0; i < m.Dim(); i++ {
		if visited[i] {
			continue
		}
		friends := make(FriendSet)
		findFriends(m, i, friends, threshold)
		visited.Union(friends)

		if len(friends) > 0 {
			group := make([]int, len(friends))
			j := 0
			for k, _ := range friends {
				group[j] = k
				j++
			}
			sort.Ints(group)
			groups = append(groups, group)
		}
	}
	return groups
}

func main() {
	var threshold float64

	flag.Float64Var(&threshold, "t", 1.0, "Threshold")
	flag.Parse()

	g := genomes.LoadGenomes(flag.Arg(0), "", false)
	groups := GroupGenomes(g, threshold)

	for i, group := range groups {
		s := make([]string, 0)
		for _, index := range group {
			s = append(s, fmt.Sprintf("%d", index))
		}
		fmt.Printf("Group %d: %s\n", i, strings.Join(s, ","))
	}
}
