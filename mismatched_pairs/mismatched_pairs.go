package main

import (
	"fmt"
	"genomics/reads"
	"genomics/utils"
	"path"
	"strings"
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

type NameSet map[string]bool

func LoadNames(fname string) NameSet {
	ret := make(NameSet)
	c := make(chan reads.ReadMsg)
	go reads.ParseFastq(fname, c)
	for {
		read := <-c
		if read.End {
			break
		}
		ret[read.Name] = true
	}
	fmt.Printf("%d SC2 matching reads\n", len(ret))
	return ret
}

func getNumber(name string) int {
	return utils.Atoi(strings.Split(name, ".")[1])
}

/*
Output the reads in fnameA that don't match well to anything in B
*/
func Unmatched(fnameA, fnameB string, onlyNames NameSet) {
	aChan := make(chan reads.ReadMsg)
	bChan := make(chan reads.ReadMsg)

	fd, fp := utils.WriteFile("unmatched.fastq")
	defer fd.Close()

	go reads.ParseFastq(fnameA, aChan)
	go reads.ParseFastq(fnameB, bChan)

	mismatched := 0
	tried := 0

	var skipA, skipB bool
	var aRead, bRead reads.ReadMsg

	for {
		if !skipA {
			aRead = <-aChan
		}
		if !skipB {
			bRead = <-bChan
		}
		skipA, skipB = false, false
		if aRead.End {
			break
		}

		if onlyNames != nil {
			if !onlyNames[aRead.Name] {
				continue
			}
		}

		if aRead.Name != bRead.Name {
			aNum := getNumber(aRead.Name)
			bNum := getNumber(bRead.Name)

			if aNum < bNum {
				skipB = true
				continue
			}
			if aNum > bNum {
				skipA = true
				continue
			}
		}

		tried++
		if !QuickAlign(aRead.Nts, utils.ReverseComplement(bRead.Nts), 1, 15) {
			mismatched++
			aRead.Output(fp)
		}
	}
	fp.Flush()
	fmt.Printf("%d/%d were mismatched pairs\n", mismatched, tried)
	fmt.Printf("Wrote unmatched.fastq")
}

func main() {
	root := "/fs/j/genomes/raw_reads/AntarcticSamples"
	fnameB := path.Join(root, "SRR13441704_2.fastq")
	fnameA := path.Join(root, "SRR13441704_1.fastq")

	// names := LoadNames(path.Join(root, "04_2-matches.fastq"))
	Unmatched(fnameA, fnameB, nil)
}
