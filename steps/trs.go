package main

import (
	"fmt"
	"genomics/utils"
)

var trsMarkers string

func init() {
	utils.Lines("./trs_locations", func(line string, err error) bool {
		location := utils.Atoi(line) - 1
		// Those are all the places where ACGAAC occurs, which is the TRS
		// marker for SARS-CoV-2 (and also SARS1, and I think all of these).
		trsMarkers = fmt.Sprintf(`%s
set arrow from %d,graph 0 to %d,graph 1 nohead filled lc "red"`,
			trsMarkers, location, location)
		return true
	})
}
