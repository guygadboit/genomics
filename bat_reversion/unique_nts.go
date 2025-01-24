package main

import (
	"fmt"
	"genomics/genomes"
)

func UniqueNts(g *genomes.Genomes) QuirkMap {
	ret := make(QuirkMap)

	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < g.Length(); j++ {
			us := g.Nts[i][j]
			unique := true
			for k := 0; k < g.NumGenomes(); k++ {
				if k == i {
					continue
				}
				if g.Nts[k][j] == us {
					unique = false
					break
				}
			}
			if unique {
				fmt.Printf("%d got %d%c\n", i, j, us)
				ret[i]++
			}
		}
	}
	return ret
}
