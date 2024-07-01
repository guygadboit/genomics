package main

import (
	"bufio"
	"fmt"
	"genomics/database"
	"genomics/genomes"
	"genomics/mutations"
	"genomics/utils"
	"log"
	"math/rand"
	"os"
	"slices"
	"strings"
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

// All the places where we found matches
type LocationSet map[int]bool

func (ls LocationSet) Print() {
	locations := utils.FromSet(ls)
	slices.Sort(locations)
	for _, pos := range locations {
		fmt.Println(pos)
	}
}

func CountOutgroupMatches(db *database.Database, ids []database.Id,
	g *genomes.Genomes, mask *genomes.Genomes, silent bool,
	description string, minOR, maxP float64) {

	odds := OutgroupOdds(g, mask)

	for _, id := range ids {
		r := &db.Records[id]
		ct, muts := CheckRecord(r, g, mask, odds, silent)
		CheckCT(ct, r, muts, description, minOR, maxP)
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

func TotalSpectrum(db *database.Database) *TransitionCounter {
	ret := NewTransitionCounter()
	from := utils.Date(2020, 1, 1)
	to := utils.Date(2020, 12, 31)

	matches := db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		coll := r.CollectionDate
		return coll.Compare(from) >= 0 && coll.Compare(to) <= 0
	})

	for match, _ := range matches {
		r := db.Records[match]
		for _, mut := range r.NucleotideChanges {
			ret.Add(mut)
		}
	}

	return ret
}

func LoadShortNames() []string {
	ret := make([]string, 0)
	utils.Lines("./short_names.txt", func(line string, err error) bool {
		fields := strings.Split(line, " ")
		ret = append(ret, fields[1])
		return true
	})
	return ret
}

func runSimulation(g *genomes.Genomes, nd *mutations.NucDistro) {
	g2 := g.Filter(0, 58)
	mask := g.Filter(5, 6, 7, 8, 10, 11)
	Simulate(g2, mask, nd, 87, 87+402, 1e4)
}

func PlotSignificant(
	db *database.Database, ids []database.Id, g *genomes.Genomes,
	minOR float64, maxP float64, silent bool,
	fname string, notes string, mask []int, expected []MatchOdds) {

	f, _ := os.Create(fname)
	defer f.Close()
	w := bufio.NewWriter(f)

	fmt.Fprintf(f, "%s. minOR=%.0f maxP=%g mask: %t silent: %t\n",
		notes, minOR, maxP, mask != nil, silent)

	individuals := CountSignificant(db, ids,
		g, minOR, maxP, silent, mask, expected)

	shortNames := LoadShortNames()
	maskSet := utils.ToSet(mask)
	for i, r := range individuals {
		if i == 0 || maskSet[i] {
			continue
		}
		name := shortNames[i]
		fmt.Fprintln(w, name, r)
	}
	w.Flush()
	fmt.Printf("Wrote %s\n", fname)
}

func main() {
	rand.Seed(9879)

	/*
		g := genomes.LoadGenomes("../fasta/more_relatives.fasta",
			"../fasta/WH1.orfs", false)
		g.RemoveGaps()
	*/
	g := genomes.LoadGenomes("./RelativesPlusKhosta.fasta",
		"../fasta/WH1.orfs", false)
	/*
		g := genomes.LoadGenomes("../fasta/all.fasta", "../fasta/WH1.orfs", false)
	*/
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

	/*
		ShowAllPossibleSilentMuts(g)
		return
	*/

	mask := []int{5, 6, 7, 8, 10, 11}
	// mask := []int{460, 461, 462, 465, 466, 463}
	// var mask []int
	db := database.NewDatabase()

	/*
		nd := mutations.NewNucDistro(g)
		runSimulation(g, nd)
		return
	*/

	cutoff := utils.Date(2020, 12, 31)
	ids := db.Filter(nil, func(r *database.Record) bool {
		if r.Host != "Human" {
			return false
		}
		/*
			if r.Country == "Egypt" {
				return false
			}
		*/
		/*
			if r.GisaidAccession != "EPI_ISL_8193636" {
				return false
			}
		*/
		if len(r.NucleotideChanges) == 0 {
			return false
		}
		return r.CollectionDate.Compare(cutoff) < 0
	})

	idSlice := utils.FromSet(ids)
	slices.Sort(idSlice)

	expected := FindAllOdds(g, mask)
	PlotSignificant(db, idSlice, g, 10, 1e-4, true,
		"individuals.txt",
		"All",
		mask, expected)
	return

	/*
		pangolins := g.Filter(0, 35, 36, 37, 38, 39, 40, 41)
		// pangolins := g.Filter(0, 33, 35)
		maskGenomes := g.Filter(mask...)
		CountOutgroupMatches(db, idSlice,
			pangolins, maskGenomes, false, "Pangolins", 10, 1e-4)
		return
	*/

	/*
		ts := TotalSpectrum(db)
		ts.Print()
		return
	*/

	/*
		fmt.Println("Before 2020-04-01")
		locations, spectrum := CountOutgroupMatches(db, pangolins, mask,
			utils.Date(2020, 1, 1),
			utils.Date(2020, 12, 31))

		locations.Print()
		spectrum.Print()
		return

		fmt.Println("After 2020-09-01")
		CountOutgroupMatches(db, nd, bats, pangolins,
			utils.Date(2020, 12, 1),
			utils.Date(2020, 12, 31))
		return
	*/

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

	// ShowOutgroupMatches(g, sc1.NucleotideChanges)
}
