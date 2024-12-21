package main

import (
	"genomics/utils"
	"fmt"
	"strings"
)

type Read struct {
	Name	string
	Score	int
}

type ReadSet struct {
	Name	string
	Reads	[]Read
}

func ParseReads(fname string) *ReadSet {
	utils.Lines(fname, func(line string, err error) bool {
		fields := strings.Split(line, "\t")
		if len(fields) < 12 {
			return true
		}
		name := fields[0]
		score := utils.Atoi(fields[11][5:])
		fmt.Println(name, score)
		return true
	})
	return nil
}

func main() {
	ParseReads("/home/benc/tmp/1")
}
