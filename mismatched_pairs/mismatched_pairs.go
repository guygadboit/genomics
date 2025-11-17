package main

import (
	"bufio"
	"genomics/reads"
	"genomics/utils"
	"os"
	"path"
)

// Can you find n in a row that match with fewer than e errors?
func QuickAlign(a, b []byte, e int, n int) bool {
	align := func(a, b []byte, n int) bool {
	offsets:
		for offset := 0; offset < len(a)-n; offset++ {
			errors := 0
			for i := 0; i < n; i++ {
				if a[i+offset] != b[i] {
					errors++
					if errors == e {
						continue offsets
					}
				}
			}
			return true
		}
		return false
	}
	return align(a, b, n) || align(b, a, n)
}

/*
Output the reads in fnameA that don't match well to anything in B
*/
func Unmatched(fnameA, fnameB string) {
	aChan := make(chan reads.ReadMsg)
	bChan := make(chan reads.ReadMsg)

	go reads.ParseFastq(fnameA, aChan)
	go reads.ParseFastq(fnameB, bChan)

	w := bufio.NewWriter(os.Stdout)

	for {
		aRead := <-aChan
		bRead := <-bChan
		if aRead.End {
			break
		}

		if !QuickAlign(aRead.Nts, utils.ReverseComplement(bRead.Nts), 1, 15) {
			aRead.Output(w)
		}
	}
	w.Flush()
}

func main() {
	root := "/fs/j/genomes/raw_reads/AntarcticSamples"
	fnameA := path.Join(root, "SRR13441704_2.fastq")
	fnameB := path.Join(root, "SRR13441704_1.fastq")
	Unmatched(fnameA, fnameB)
}
