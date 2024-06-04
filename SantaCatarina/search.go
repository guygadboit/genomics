package main

import (
	"bufio"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/stats"
	"genomics/utils"
	"log"
	"os"
	"time"
)

func SearchDB(db *database.Database) {
	// interesting := database.ParseMutations("T5929G,A8651C,G16206A")
	// interesting := database.ParseMutations("T5929G,T8601C,A8651C,G16206A,T19218G")

	all := database.ParseMutations(
		"G3392C,A5706G,T5929G,G7587A,A8576C,T8601C,A8651C,C14408T,C14422T," +
			"T14895G,A15106C,T15132A,G16206A,T16251C,T19218G,C19273G," +
			"G22773A,A23403G")

	matches := make([]database.Id, 0)
	results := make(map[database.Id]database.Mutations)

	for _, r := range db.Records {
		got := r.HasMuts(all)
		if len(got) != 0 {
			matches = append(matches, r.Id)
			results[r.Id] = got
		}
	}
	db.Sort(matches, database.COLLECTION_DATE)

	for _, id := range matches {
		r := db.Records[id]
		got := results[r.Id]
		fmt.Printf("%d %s %s %s %s %s %s %s %s %d\n", len(got),
			r.GisaidAccession, r.CollectionDate.Format(time.DateOnly),
			r.SubmissionDate.Format(time.DateOnly),
			r.Country, r.Region, r.City, r.Host,
			got.ToString(), len(r.NucleotideChanges))
	}
}

func showEarly(db *database.Database) {
	var v Visualizer
	v.Init()
	cutoff := utils.Date(2021, 1, 1)

	v.Special("C14408T", "*")
	v.Special("A23403G", "#")
	v.Special("A5706G", "@")

	v.Special("T5929G", "p1")
	v.Special("A8651G", "p2")
	v.Special("G16206A", "p3")

	v.Special("G3392C", "o1")
	v.Special("G7587A", "o2")
	v.Special("A8576C", "o3")
	v.Special("T8601C", "o4")
	v.Special("C14422T", "o5")
	v.Special("T14895G", "o6")
	v.Special("A15106C", "o7")
	v.Special("T5132A", "o8")
	v.Special("T6251C", "o9")
	v.Special("T9218G", "o10")
	v.Special("C9273G", "o11")
	v.Special("G22773A", "o12")

	for _, r := range db.Records {
		if r.CollectionDate.Compare(cutoff) < 0 {
			if r.Host == "Human" {
				muts := string(v.Encode(r.NucleotideChanges))
				fmt.Printf("%s %s %s %s %s %s\n", r.GisaidAccession,
					r.CollectionDate.Format(time.DateOnly),
					r.Country, r.Region, r.City, muts)
			}
		}
	}
	v.Show()
}

func CountOutgroupMatches(db *database.Database, nd *mutations.NucDistro,
	bats *genomes.Genomes,
	pangolins *genomes.Genomes,
	from time.Time, to time.Time) {

	its := 1000
	batHits := OutgroupMontecarlo(bats, nd, its)
	pangHits := OutgroupMontecarlo(pangolins, nd, its)

	matches := db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		coll := r.CollectionDate
		return coll.Compare(from) >= 0 && coll.Compare(to) <= 0
	})
	fmt.Printf("Considering %d sequences\n", len(matches))

	ids := utils.FromSet(matches)
	db.Sort(ids, database.COLLECTION_DATE)

	checkMatches := func(r *database.Record,
		matches Matches, hits int, tag string) {
		if len(matches) == 0 {
			return
		}
		var ct stats.ContingencyTable
		n := len(matches)
		ct.Init(n, len(r.NucleotideChanges)-n, hits, its-hits)
		OR, p := ct.FisherExact()

		if p < 1e-3 {
			fmt.Printf("%s %d/%d match %s: %s OR=%f p=%.4g\n", r.ToString(),
				n, len(r.NucleotideChanges), tag, matches.ToString(false),
				OR, p)
		}
	}

	for _, id := range ids {
		r := &db.Records[id]

		batMatches := OutgroupMatches(bats,
			r.NucleotideChanges, "b", false)
		pangMatches := OutgroupMatches(pangolins,
			r.NucleotideChanges, "p", false)

		batMatches, pangMatches = RemoveIntersection(batMatches, pangMatches)

		checkMatches(r, batMatches, batHits, "Bat")
		checkMatches(r, pangMatches, pangHits, "Pangolin")
	}
}

func mutsDates(db *database.Database) {
	f, _ := os.Create("muts-dates.txt")
	defer f.Close()
	w := bufio.NewWriter(f)

	for _, r := range db.Records {
		d := r.CollectionDate.Unix()
		if r.Host == "Human" && d > 0 {
			numMuts := len(r.NucleotideChanges)
			fmt.Fprintf(w, "%d %s %s %d\n", d, r.Country,
				r.GisaidAccession, numMuts)
		}
	}

	w.Flush()
	fmt.Println("Wrote muts-dates.txt")
}

func main() {
	g := genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
		"../fasta/WH1.orfs", false)
	/*
		// OutgroupMontcarlo(g, 1000, 17)
		ShowOutgroupMatches(g)
	*/

	/*
		var db database.Database
		db.Parse(database.ROOT + "gisaid2020.tsv.gz")
		db.Save(database.GOB_NAME)
		return
	*/

	db := database.NewDatabase()
	nd := mutations.NewNucDistro(g)

	pangolins := g.Filter(0, 35, 36, 37, 38, 39, 40, 41)
	bats := g.Filter(0, 8, 5, 6, 7, 10, 11)

	fmt.Println("Before 2020-03-01")
	CountOutgroupMatches(db, nd, bats, pangolins,
		utils.Date(2020, 1, 1),
		utils.Date(2020, 3, 1))

	fmt.Println("After 2020-09-01")
	CountOutgroupMatches(db, nd, bats, pangolins,
		utils.Date(2020, 12, 1),
		utils.Date(2020, 12, 31))
	return

	/*
		mutsDates(db)
		return
	*/

	muts := database.ParseMutations("A5706G")
	br := db.SearchByMuts(muts, 1)
	br = db.Filter(br, func(r *database.Record) bool {
		return r.CollectionDate.Year() == 2020 && r.Country == "Brazil"
	})

	// Now let's add anything early than EPI_ISL_541370, which is that first
	// Human Santa Catarina sample, which has some of the same muts, regardless
	// of location.
	sc1 := &db.Records[db.GetByAccession("EPI_ISL_541370")[0]]
	others := db.SearchByMuts(sc1.NucleotideChanges, 7)
	others = db.Filter(others, func(r *database.Record) bool {
		cutoff := utils.Date(2020, 3, 1)
		return r.CollectionDate.Compare(cutoff) <= 0
	})

	// And now anything very early, regardless of muts
	early := db.Filter(nil, func(r *database.Record) bool {
		cutoff := utils.Date(2020, 1, 10)
		return r.CollectionDate.Compare(cutoff) <= 0 && r.Host == "Human"
	})

	utils.Union(others, early)
	records := utils.FromSet(others)
	db.Sort(records, database.COLLECTION_DATE)

	// Now put the BR ones at the start
	brRecords := utils.FromSet(br)
	db.Sort(brRecords, database.COLLECTION_DATE)
	records = append(brRecords, records...)

	f, err := os.Create("table.html")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	MakeTable(db, g, records, sc1.NucleotideChanges, w)
	w.Flush()
	fmt.Printf("Wrote table.html\n")

	ShowOutgroupMatches(g, sc1.NucleotideChanges)
}
