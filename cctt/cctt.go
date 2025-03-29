package main

import (
	"flag"
	"fmt"
	"genomics/database"
	"genomics/pileup"
	"genomics/utils"
	"log"
	"os"
	"path"
	"time"
)

type Contents struct {
	C8782 int // The depth of C at 8782
	T8782 int

	C28144 int
	T28144 int

	Classification string
}

func (c *Contents) ToString() string {
	return fmt.Sprintf("%s %d|%d %d|%d", c.Classification,
		c.C8782, c.T8782, c.C28144, c.T28144)
}

func Min(x, y int) int {
	if x < y {
		return x
	} else {
		return y
	}
}

func Classify(pu *pileup.Pileup, minDepth int) Contents {
	pos8782 := pu.Get(8782)
	pos28144 := pu.Get(28144)

	ret := Contents{
		pos8782.GetDepthOf('C'), pos8782.GetDepthOf('T'),
		pos28144.GetDepthOf('C'), pos28144.GetDepthOf('T'),
		"*",
	}

	if ret.C8782 > minDepth &&
		ret.T8782 < minDepth &&
		ret.C28144 > minDepth &&
		ret.T28144 < minDepth {
		ret.Classification = "CC"
	} else if ret.T8782 > minDepth &&
		ret.C8782 < minDepth &&
		ret.T28144 > minDepth &&
		ret.C28144 < minDepth {
		ret.Classification = "TT"
	}

	return ret
}

func FindPileup(record *database.Record, root string) *pileup.Pileup {
	path := path.Join(root, fmt.Sprintf("%s-WH1-index.txt.gz", record.SRAs()))
	if _, err := os.Stat(path); err != nil {
		log.Printf("%s doesn't exist\n", path)
		return nil
	}
	ret, err := pileup.Parse2(path)
	if err != nil {
		log.Fatal(err)
	}
	return ret
}

func LoadRecords(db *database.Database, fname string) []database.Id {
	ret := make([]database.Id, 0)
	utils.Lines(fname, func(line string, lineErr error) bool {
		ids := db.GetByAccession(line)
		ret = append(ret, ids...)
		return true
	})
	return ret
}

func AnalyseReads(db *database.Database,
	ids []database.Id, minDepth int, prefix string) {
	fmt.Println("Date AccNo SRA Country class 8782:C|T 28144:C|T")

	root := "/fs/bowser/genomes/raw_reads/"
	for _, id := range ids {
		record := db.Get(id)
		pu := FindPileup(record, path.Join(root, prefix))
		if pu == nil {
			continue
		}
		contents := Classify(pu, minDepth)
		fmt.Println(record.CollectionDate.Format(time.DateOnly),
			record.GisaidAccession, record.SRAs(),
			record.Country, contents.ToString())
	}
}

func main() {
	var minDepth int

	flag.IntVar(&minDepth, "min-depth", 3, "Min depth")
	flag.Parse()

	db := database.NewDatabase()

	cc := LoadRecords(db, "./cc")
	AnalyseReads(db, cc, minDepth, "CC")

	/*
		tt := LoadRecords(db, "./tt")
		AnalyseReads(db, tt, minDepth, "TT")
	*/
}
