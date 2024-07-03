package main

import (
	"strings"
	"fmt"
	"genomics/utils"
	"genomics/genomes"
	"bufio"
	"log"
	"os"
)

type Read struct {
	Name	string
	Guff	string
	Nts     []byte
	Quality []byte
}

type ReadMsg struct {
	Read
	End	bool
}

func (r *Read) Output(w *bufio.Writer) {
	fmt.Fprintf(w, "@%s\n", r.Guff)
	fmt.Fprintf(w, "%s\n", string(r.Nts))
	fmt.Fprintf(w, "+%s\n", r.Guff)
	fmt.Fprintf(w, "%s\n", string(r.Quality))
}

func ParseFastq(fname string, output chan ReadMsg) error {
	var msg ReadMsg
	var destination *[]byte
	var post bool
	var ret error

	utils.Lines(fname, func(line string, err error) bool {
		if err != nil {
			ret = err
			return false
		}
		if line[0] == '@' {
			fields := strings.Split(line, " ")
			msg.Read.Name = fields[0][1:]
			msg.Read.Guff = line[1:]
			destination = &msg.Read.Nts
		} else if line[0] == '+' {
			destination = &msg.Read.Quality
			post = true
		} else {
			*destination = []byte(line)
			if post {
				output <- msg
				post = false
			}
		}
		return true
	})
	if ret != nil {
		return ret
	}
	output <- ReadMsg{Read{"", "", nil, nil}, true}
	return nil
}

func test() {
	// nts := []byte("GCTCAAAGGAGTCAAATTACATTACACATAAACGAACTTATGGATTTGTTTATG")

	nts := []byte("GTGCTTGCATACGTAGACCATTCTTATGTTGTAAATGCTGTTACGACCATGTCATATCAACATCACATAAATTAGTCTTGTCTGTTAATCCGTATGTTTGCAATGCTCCAGGTTGTGATGTCACAGATGTGACTCAACTTTACTTAGGAGG")

	// wh1 := "/fs/f/genomes/viruses/SARS2/index"
	hn2021 := "/fs/f/genomes/viruses/HN2021/HN2021A-index"

	var search genomes.IndexSearch
	var inWH1 bool

	for search.Init(hn2021, nts); !search.End(); search.Next() {
		inWH1 = true
		break
	}
	fmt.Println(inWH1)
}

func main() {
	reads := make(chan ReadMsg)
	go ParseFastq("SRR29285436.fastq.gz", reads)

	wh1 := "/fs/f/genomes/viruses/SARS2/index"

	var notFoundCount, total int

	f, err := os.Create("interesting.fastq")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)

	for {
		read := <- reads
		if read.End {
			break
		}

		var search genomes.BidiIndexSearch
		var found bool

		for search.Init(wh1, read.Nts); !search.End(); search.Next() {
			found = true
			break
		}
		
		if !found {
			notFoundCount++
			read.Output(w)
		}
		total++

		if total % 10000 == 0 {
			fmt.Println(total)
		}
	}
	fmt.Printf("Not perfect matches in WH1: %d/%d\n", notFoundCount, total)
	w.Flush()
}
