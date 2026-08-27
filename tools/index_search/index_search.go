package main

import (
	"flag"
	"fmt"
	"genomics/genomes"
	"log"
)

func Search(root string, needle []byte) {
	count := 0
	for search := genomes.NewBidiIndexSearch(root,
		needle); !search.End(); search.Next() {
		pos, err := search.Get()
		if err != nil {
			log.Fatal(err)
		}
		var f string
		if search.IsForwards() {
			f = "+"
		}
		fmt.Printf("%d%s\n", pos+1, f)
		count++
	}
	fmt.Printf("%d matches\n", count)
}

func main() {
	var pattern string
	flag.StringVar(&pattern, "p", "", "Pattern")
	flag.Parse()

	needle := []byte(pattern)
	for _, arg := range flag.Args() {
		fmt.Println(arg)
		Search(arg, needle)
	}
}
