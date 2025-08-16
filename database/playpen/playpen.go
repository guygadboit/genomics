package main

import (
	"fmt"
	"genomics/database"
	"genomics/stats"
	"genomics/utils"
	"slices"
	"strings"
	"time"
)

var DiamondPrincess []string = []string{
	"EPI_ISL_416591",
	"EPI_ISL_416633",
	"EPI_ISL_416593",
	"EPI_ISL_416585",
	"EPI_ISL_416570",
	"EPI_ISL_416578",
	"EPI_ISL_416600",
	"EPI_ISL_416575",
	"EPI_ISL_416584",
	"EPI_ISL_416631",
	"EPI_ISL_416576",
	"EPI_ISL_412968",
	"EPI_ISL_416614",
	"EPI_ISL_416581",
	"EPI_ISL_416625",
	"EPI_ISL_416573",
	"EPI_ISL_416583",
	"EPI_ISL_416603",
	"EPI_ISL_416571",
	"EPI_ISL_416620",
	"EPI_ISL_416588",
	"EPI_ISL_416589",
	"EPI_ISL_416613",
	"EPI_ISL_416577",
	"EPI_ISL_416604",
	"EPI_ISL_416612",
	"EPI_ISL_416624",
	"EPI_ISL_416608",
	"EPI_ISL_416610",
	"EPI_ISL_416596",
	"EPI_ISL_416579",
	"EPI_ISL_416572",
	"EPI_ISL_416601",
	"EPI_ISL_416634",
	"EPI_ISL_416580",
	"EPI_ISL_454749",
	"EPI_ISL_416605",
	"EPI_ISL_416621",
	"EPI_ISL_416574",
	"EPI_ISL_416569",
	"EPI_ISL_416565",
	"EPI_ISL_416630",
	"EPI_ISL_416594",
	"EPI_ISL_416609",
	"EPI_ISL_416597",
	"EPI_ISL_416628",
	"EPI_ISL_412969",
	"EPI_ISL_416611",
	"EPI_ISL_416607",
	"EPI_ISL_416582",
	"EPI_ISL_416567",
	"EPI_ISL_416566",
	"EPI_ISL_416587",
	"EPI_ISL_416599",
	"EPI_ISL_416595",
	"EPI_ISL_416606",
	"EPI_ISL_416590",
	"EPI_ISL_416586",
	"EPI_ISL_420889",
	"EPI_ISL_416619",
	"EPI_ISL_416622",
	"EPI_ISL_416592",
	"EPI_ISL_416617",
	"EPI_ISL_416598",
	"EPI_ISL_416627",
	"EPI_ISL_416629",
	"EPI_ISL_416632",
}

var Outliers []string = []string{
	"EPI_ISL_406592",
	"EPI_ISL_414588",
	"EPI_ISL_412900",
	"EPI_ISL_408487",
	"EPI_ISL_408483 ",
	"EPI_ISL_406595",
}

func interestingMuts(r *database.Record) string {
	muts := database.ParseMutations("T18060C,T8782C,C28144T," +
		"C28657T,T9477A,C28863T,G25979T")

	interesting := map[database.Mutation]string{
		muts[0]: "α1",
		muts[1]: "α2",
		muts[2]: "α3",
		muts[3]: "α1a",
		muts[4]: "α1b",
		muts[5]: "α1c",
		muts[6]: "α1d",
	}
	ret := make([]string, 0)
	for _, mut := range r.NucleotideChanges {
		mut.Silence = utils.UNKNOWN
		mut.From, mut.To = mut.To, mut.From
		name, there := interesting[mut]
		if there {
			ret = append(ret, name)
		}
	}
	return strings.Join(ret, ",")
}

type RdRPMutation struct {
	database.Mutation
	ids []database.Id
}

func RdRPVariants(db *database.Database) {
	// The RdRP is nsp7, nsp8 and nsp12. Find all sequences with variations
	// anywhere in there.
	positions := make([]utils.OneBasedPos, 0)

	addRange := func(start, end utils.OneBasedPos) {
		for i := start; i <= end; i++ {
			positions = append(positions, i)
		}
	}

	addRange(11843, 12091) // nsp7
	addRange(12094, 12685) // nsp8
	addRange(13442, 16236) // nsp12

	posSet := utils.ToSet(positions)
	matches := db.SearchByMutPosition(positions, 1)

	// Now group them by the actual unique mutations. Make an index of mutation
	// to ids that have it.
	index := make(map[database.Mutation]database.IdSet)

	for _, match := range matches {
		r := db.Get(match.Id)
		for _, mut := range r.NucleotideChanges {

			// Don't print out anything referring to other muts these records
			// may have that aren't in nsp7/8/12
			if !posSet[mut.Pos] {
				continue
			}
			_, there := index[mut]
			if !there {
				index[mut] = make(database.IdSet)
			}
			index[mut][match.Id] = true
		}
	}

	// Now make that into something we can sort
	muts := make([]RdRPMutation, 0, len(index))
	for k, v := range index {
		muts = append(muts, RdRPMutation{k, utils.FromSet(v)})
	}

	slices.SortFunc(muts, func(a, b RdRPMutation) int {
		ka, kb := len(a.ids), len(b.ids)
		if ka < kb {
			return 1
		}
		if ka > kb {
			return -1
		}
		return 0
	})

	for _, mut := range muts {
		fmt.Printf("%s in %d sequences\n", mut.ToString(), len(mut.ids))

		db.Sort(mut.ids, database.COLLECTION_DATE)
		for _, id := range mut.ids {
			r := db.Get(id)
			fmt.Println(r.ToString())
		}
	}
}

func Pangolin(db *database.Database) {
	positions := []utils.OneBasedPos{8782}
	matches := db.SearchByMutPosition(positions, 1)
	for _, m := range matches {
		var matched bool
		r := db.Get(m.Id)
		if r.Host != "Human" {
			continue
		}
		for _, mut := range r.NucleotideChanges {
			if mut.From == 'C' && mut.To == 'T' {
				matched = true
				break
			}
		}
		if !matched {
			continue
		}
		fmt.Println(r.ToString())
	}
}

func D614(db *database.Database) {
	ids := utils.FromSet(db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		for _, m := range r.AAChanges {
			if m.Gene == "S" && m.Pos == 614 && m.From == 'D' && m.To == 'G' {
				return false
			}
		}
		return true
	}))

	db.Sort(ids, database.COLLECTION_DATE)

	for _, id := range ids {
		r := db.Get(id)
		fmt.Printf("%s %d\n", r.ToString(), len(r.AAChanges))
	}
}

func TT(db *database.Database) {
	ids := utils.FromSet(db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}

		if r.CollectionDate.Compare(utils.Date(2020, 2, 1)) < 0 {
			return false
		}

		/*
			if len(r.NucleotideChanges) != 0 {
				return false
			}
		*/

		/*
			for _, m := range r.AAChanges {
				if m.Gene == "S" && m.Pos == 614 && m.From == 'D' && m.To == 'G' {
					return false
				}
			}
		*/

		good := false
		for _, m := range r.NucleotideChanges {
			// if m.Pos == 8782 && m.To == 'T' {
			if m.Pos == 8878 && m.To == 'T' {
				good = true
			}
			if m.Pos == 28144 {
				good = false
			}
		}
		return good
	}))

	db.Sort(ids, database.COLLECTION_DATE)

	for _, id := range ids {
		r := db.Get(id)
		fmt.Printf("%s\n", r.Summary())
	}
	fmt.Printf("%d sequences with the required filters\n", len(ids))
}

func CTRate(db *database.Database) {
	ct, nonCt := 0, 0

	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}

		if r.CollectionDate.Compare(utils.Date(2020, 2, 1)) < 0 {
			return false
		}

		for _, m := range r.NucleotideChanges {
			if m.From == 'C' && m.To == 'T' {
				ct++
			} else {
				nonCt++
			}
		}
		return false
	})

	fmt.Printf("CT: %d nonCt: %d\n", ct, nonCt)
}

func NoMuts(db *database.Database) {
	regions := make(map[string]*stats.ContingencyTable)
	C, D := 0, 0

	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		/*
			if r.CollectionDate.Compare(utils.Date(2020, 6, 1)) < 0 {
				return false
			}
		*/

		var matched bool
		// matched = len(r.NucleotideChanges) == 0
		matched = len(r.NucleotideChanges) == 1

		/*
			if len(r.NucleotideChanges) == 1 {
				ref := database.Mutation{23403, 'A', 'G', utils.NON_SILENT}
				if r.NucleotideChanges[0] == ref {
					matched = true
				}
			}
		*/

		region := r.Region
		ct, there := regions[region]
		if !there {
			var newCt stats.ContingencyTable
			regions[region] = &newCt
			ct = regions[region]
		}
		if matched {
			C++
			ct.A++
		} else {
			ct.B++
			D++
		}

		if matched && len(r.SRA) > 1 {
			fmt.Println(r.Summary())
		}

		return false
	})

	type result struct {
		region string
		ct     stats.ContingencyTable
	}
	results := make([]result, 0)

	for k, _ := range regions {
		regions[k].C = C
		regions[k].D = D
		if regions[k].A > 0 {
			ct := *regions[k]
			ct.FisherExact(stats.GREATER)
			results = append(results, result{k, ct})
		}
	}

	slices.SortFunc(results, func(a, b result) int {
		if a.ct.P < b.ct.P {
			return -1
		}
		if a.ct.P < b.ct.P {
			return 1
		}
		return 0
	})

	for _, r := range results {
		fmt.Println(r.region, r.ct.A, r.ct.B, r.ct.OR, r.ct.P)
	}

}

func EarlyLineages(db *database.Database) {
	counts := make(map[string]int)
	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}

		if len(r.SRA) == 0 {
			return false
		}

		/*
			if r.CollectionDate.Compare(utils.Date(2020, 7, 1)) > 0 {
				return false
			}
		*/

		/*
			if len(r.NucleotideChanges) > 2 {
				return false
			}
		*/

		C1, C2 := true, false
		allowed := 0
		for _, m := range r.NucleotideChanges {
			if m.Pos == 8782 && m.To == 'T' {
				C1 = false
				allowed++
			}
			if m.Pos == 28144 && m.To == 'C' {
				C2 = true
				allowed++
			}
		}

		/*
			if len(r.NucleotideChanges) != allowed {
				return false
			}
		*/

		/*
			if len(r.NucleotideChanges) < allowed+5 {
				return false
			}
		*/

		var class string

		display := func() {
			fmt.Println(class, r.SRAs(), len(r.NucleotideChanges),
				r.GisaidAccession, r.Country,
				r.CollectionDate.Format(time.DateOnly))
		}

		if C1 {
			if C2 {
				class = "CC"
				display()
			} else {
				class = "CT"
				// fmt.Println(r.GisaidAccession)
			}
		} else {
			if C2 {
				class = "TC"
				// fmt.Println(r.GisaidAccession)
			} else {
				class = "TT"
				display()
			}
		}
		counts[class]++
		// fmt.Println(class, r.Summary())
		return false
	})

	fmt.Println(counts)
}

func EarlyReads(db *database.Database) {
	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		if len(r.SRA) == 0 {
			return false
		}
		sras := strings.Join(r.SRA, ",")

		fmt.Println(sras, len(r.NucleotideChanges),
			r.GisaidAccession, r.SampleSRA, r.Country,
			r.CollectionDate.Format(time.DateOnly))
		return true
	})
}

func NonHuman(db *database.Database) {
	db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" && len(r.SRA) != 0 {
			fmt.Println(r.SRAs(), r.GisaidAccession, r.Host)
		}
		return true
	})
}

func ORF8LI(db *database.Database) {
	keys := []database.AAMutationPos{database.AAMutationPos{"ORF8", 10}}

	s := database.NewAAMutationIndexSearch(db, keys)
	results, _ := s.Get(0)
	for id, _ := range results {
		record := db.Get(id)
		if record.Host != "Human" {
			continue
		}
		for _, mut := range record.AAChanges {
			if mut.Pos == 10 && mut.To == 'L' {
				fmt.Println(db.Get(id).Summary())
			}
		}
	}
}

func main() {
	db := database.NewDatabase()
	// EarlyReads(db)
	// NonHuman(db)
	// EarlyLineages(db)
	// NoMuts(db)
	// RdRPVariants(db)
	// Pangolin(db)
	// TT(db)
	// CTRate(db)
	// CTDistribution(db)
	ORF8LI(db)
}
