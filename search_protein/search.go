package main

import (
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"io"
	"log"
	"strings"
)

/*
	Look for len long substrings of needle in haystack
*/
func search(haystack *genomes.Genomes, needle []byte, length int) {
	for start := 0; start < len(needle); start++ {
		end := start + length
		if end > len(needle) {
			break
		}

		for i := 0; i < haystack.NumGenomes(); i++ {
			var s genomes.Search
			for s.Init(haystack, i,
				needle[start:end], 0.0); !s.End(); s.Next() {
				pos, _ := s.Get()
				fmt.Printf("Found %s at %d in %s (length %d)\n",
					string(needle[start:end]), pos, haystack.Names[i], length)
			}
		}
	}
}

func main() {
	gs := genomes.LoadGenomes(
		"/fs/f/genomes/human/GRCh38_latest_protein.faa.gz", "", false)

	for length := 7; length < 12; length++ {
		f := utils.NewFileReader("peptides.txt")
		defer f.Close()

	loop:
		for {
			line, err := f.ReadString('\n')
			switch err {
			case io.EOF:
				break loop
			case nil:
				break
			default:
				log.Fatal("Can't read peptides file")
			}
			line = strings.TrimSpace(line)
			fmt.Printf("Searching %s length %d\n", line, length)
			search(gs, []byte(line), length)
		}
	}
}
