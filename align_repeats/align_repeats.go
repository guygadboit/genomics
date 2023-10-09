package main

import (
	"bufio"
	"flag"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"os"
)

func main() {
	var (
		fasta, index, pattern string
		output                string
		context               int
		infer                 bool
	)

	// If you're using multiple names (positional args) don't use -f or -i
	flag.StringVar(&fasta, "f", "", "Fasta File (or inferred from name)")
	flag.StringVar(&index, "i", "", "Index Directory (or inferred from name)")
	flag.StringVar(&pattern, "p", "", "Pattern")
	flag.IntVar(&context, "c", 0, "How much context")
	flag.BoolVar(&infer, "infer", true, "Infer paths to fasta/index")
	flag.StringVar(&output, "o", "output.fasta", "Output file name")
	flag.Parse()

	fd, err := os.Create(output)
	if err != nil {
		log.Fatal("Can't create output file.")
	}
	defer fd.Close()
	w := bufio.NewWriter(fd)

	args := flag.Args()
	for i := 0; i < len(args); i++ {
		name := args[i]

		if infer {
			fasta = fmt.Sprintf("../fasta/bacteria/%s/%s.fasta.gz", name, name)
			index = fmt.Sprintf("../fasta/bacteria/%s/index", name)
		}

		needle := []byte(pattern)
		g := genomes.LoadGenomes(fasta, "", true)
		nts := g.Nts[0]

		var search genomes.IndexSearch
		var found bool

		for search.Init(index, needle); !search.End(); search.Next() {
			found = true
			pos, _ := search.Get()

			start := pos - context
			if start < 0 {
				start = 0
			}
			end := pos + context
			if end > g.Length() {
				end = g.Length()
			}

			fmt.Fprintf(w, ">%s-%d-%d\n", name, start, end)
			utils.Wrap(w, nts[start:end])
		}
		if found {
			fmt.Printf("Yes found in %s\n", name)
		} else {
			fmt.Printf("Not found in %s\n", name)
		}
	}

	w.Flush()
	fmt.Printf("Wrote %s\n", output)
}
