package main

import (
    "fmt"
    "genomics/genomes"
    "genomics/comparison"
)

func main() {
    g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
        "../fasta/WH1.orfs", false)

    for i := 0; i < g.NumGenomes(); i++ {
        for j := i+1; j < g.NumGenomes(); j++ {
            c := comparison.Compare(g, i, j)
            S, NS, _ := c.SilentCount()
            ratio := float64(NS)/float64(S)
            fmt.Printf("%.2f %d %d %d-%d\n", ratio, NS, S, i, j)
        }
    }
}
