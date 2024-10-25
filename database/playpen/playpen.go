package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
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
	muts := database.ParseMutations("T18060C,T8782C,C28144T,"+
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
		mut.Silence = database.UNKNOWN
		mut.From, mut.To = mut.To, mut.From
		name, there := interesting[mut]
		if there {
			ret = append(ret, name)
		}
	}
	return strings.Join(ret, ",")
}

func main() {
	db := database.NewDatabase()

	// dp := utils.ToSet(DiamondPrincess)
	muts := database.ParseMutations("T18060C,T8782C,C28144T,"+
		"C28657T,T9477A,C28863T,G25979T")
	// muts := database.ParseMutations("A17858G,C17747T")	// ν1 and ν2
	// ids := db.SearchByMuts(muts, len(muts))
	ids := db.SearchByMuts(muts, 1)
	fmt.Printf("Found %d with those muts\n", len(ids))

	// cutoff := utils.Date(2020, 1, 30)
	// outliers := utils.ToSet(Outliers)
	/*
	ids = db.Filter(ids, func(r *database.Record) bool {
		return r.Host == "Human" // && outliers[r.GisaidAccession]
		// return r.CollectionDate.Compare(cutoff) < 0
	})
	*/

	sorted := utils.FromSet(ids)
	db.Sort(sorted, database.COLLECTION_DATE)

	for _, id := range sorted {
		r := &db.Records[id]
		// fmt.Println(r.Summary(), interestingMuts(r))
		fmt.Println(r.CollectionDate.Format(time.DateOnly),
			r.SubmissionDate.Format(time.DateOnly),
			r.GisaidAccession, r.Country, r.Region,
			r.DeletionsSummary(), r.Host, interestingMuts(r))
	}
}
