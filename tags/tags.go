package main

import (
	"bufio"
	"encoding/gob"
	"fmt"
	"genomics/genomes"
	"genomics/utils"
	"log"
	"math/rand"
	"os"
	"slices"
    "reflect"
)

// A short pattern of nts, probably 6, with a lot of silent muts in it, that
// might be useful as a classifer
type Pattern struct {
	Origin int // the genome the nts we record here were in
	Pos    int
	Nts    []byte
}

// Represents the "distance" (number of silent muts) from a particular pattern
type TagDist struct {
	// The pattern that we have and the one that we are numMuts from
	Us, Them string
	NumMuts  int
}

// Patterns indexed by their nts converted to string
type PatternSet map[string]Pattern

type Tag struct {
	genomes *genomes.Genomes
	Pos     int

	// The set of alternative patterns found in all the genomes for this tag.
	Patterns PatternSet

	// Maps each genome to its distance from the tag it's closest to
	Mapping map[int]TagDist
}

// Return a map of the alleles and their counts
func (t *Tag) Alleles() map[string]int {
	ret := make(map[string]int)

	for _, v := range t.Mapping {
		ret[v.Us]++
	}

	return ret
}

func (t *Tag) Length() int {
	for k, _ := range t.Patterns {
		return len(k)
	}
	log.Fatal("Tag has no patterns!")
	return 0
}

func (t *Tag) Print() {
	fmt.Printf("Tag at %d (%s)\n", t.Pos, OrfNames.GetName(t.Pos))

	for k, v := range t.Mapping {
		fmt.Printf("%d(%s): %d has %s (%d muts from %s)\n", t.Pos,
			OrfNames.GetName(t.Pos), k, v.Us, v.NumMuts, v.Them)
	}
}

func (t *Tag) NumGroups() int {
	alleles := t.Alleles()
	return len(alleles)
}

// Does tag contain the specified genome?
func (t *Tag) Contains(which int) bool {
	for k, _ := range t.Mapping {
		if k == which {
			return true
		}
	}
	return false
}

// Show which other genomes have the same allele, or differ by only one, from
// which, sorted by index.
func (t *Tag) ShowClosest(which int) {
	friends := make([]int, 0)

	us := []byte(t.Mapping[which].Us)

	for k, td := range t.Mapping {
		nm := utils.NumMuts(us, []byte(td.Us))
		if nm <= 1 {
			friends = append(friends, k)
		}
	}
	slices.Sort(friends)

	fmt.Printf("%d: %d (%s): ", which, t.Pos, OrfNames.GetName(t.Pos))
	for _, f := range friends {
		fmt.Printf("%d ", f)
	}
	fmt.Printf("\n")
}

// Show what the specified genome has at the location where t is.
func (t *Tag) Compare(which int) {
	t.Print()
	n := t.Length()
	fmt.Printf("%d has %s\n",
		which, string(t.genomes.Nts[which][t.Pos:t.Pos+n]))
}

/*
Find Patterns that could potentially be used as tags because they are
length long and differ between two genomes by at least minMuts. Another
function will build a tag from a pattern.
*/
func FindPatterns(g *genomes.Genomes, length int, minMuts int) []Pattern {
	ret := make([]Pattern, 0)

positions:
	for i := 0; i < g.Length()-length; i++ {
		for j := 0; j < g.NumGenomes(); j++ {
			for k := 0; k < j; k++ {
				err, silent, numMuts := genomes.IsSilent(g, i, length, j, k)
				if err != nil {
					continue
				}
				if silent && numMuts >= minMuts {
					// fmt.Printf("Found tag at %d\n", i)
					ret = append(ret, Pattern{j, i, g.Nts[j][i : i+length]})
					continue positions
				}
			}
		}
	}
	return ret
}

func (t *Tag) Init(g *genomes.Genomes, p Pattern) {
	n := len(p.Nts)

	t.genomes = g
	t.Pos = p.Pos
	t.Patterns = make(PatternSet)
	t.Mapping = make(map[int]TagDist)

	// The first alternative at this location is the one we started with-- p,
	// which has a zero distance to itself.
	sp := string(p.Nts)
	t.Patterns[sp] = p
	t.Mapping[p.Origin] = TagDist{sp, sp, 0}

	for i := 0; i < g.NumGenomes(); i++ {
		// There shouldn't be an error here because we wouldn't be here if this
		// location hadn't been OK in FindPatterns
		_, silent, numMuts := genomes.IsSilent(g, p.Pos, n, p.Origin, i)

		if silent {
			alt := g.Nts[i][p.Pos : p.Pos+n]
			us := string(alt)
			t.Patterns[us] = Pattern{i, p.Pos, alt}
			t.Mapping[i] = TagDist{us, sp, numMuts}
		}
	}
}

func CreateTags(g *genomes.Genomes, patterns []Pattern) []Tag {
	ret := make([]Tag, len(patterns))
	for i, p := range patterns {
		ret[i].Init(g, p)
	}
	return ret
}

func SaveTags(fname string, tags []Tag) {
	fd, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	enc := gob.NewEncoder(fp)
	err = enc.Encode(tags)
	if err != nil {
		log.Fatal(err)
	}
}

/*
Load serialized tags and associate them with a set of genomes (which should be
identical to those used when they were originally created)
*/
func LoadTags(g *genomes.Genomes, fname string) []Tag {
	ret := make([]Tag, 0)

	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)
	dec := gob.NewDecoder(fp)
	err = dec.Decode(&ret)
	if err != nil {
		log.Fatal(err)
	}

	for i, _ := range ret {
		ret[i].genomes = g
	}

	return ret
}

// Show which triples "which" is in, and return the triples
func TriplesWith(tags []Tag, which int) {
	for _, tag := range tags {
		if tag.Contains(which) && tag.NumGroups() <= 3 {
			tag.ShowClosest(which)
		}
	}
}

// Just find the tag at a particular position
func FindTag(tags []Tag, pos int) *Tag {
	for _, t := range tags {
		if t.Pos == pos {
			return &t
		}
	}
	return nil
}

// Find tags where only one genome has a particular allele
func FindUnique(tags []Tag) {
	// Count of how many unique tags each genome has
	unique := make(map[int]int)

	for _, t := range tags {
		alleles := t.Alleles()
		for k, v := range alleles {
			if v != 1 {
				continue
			}
			for genome, td := range t.Mapping {
				if td.Us == k {
					fmt.Printf("%s at %d(%s) unique in %d\n",
						k, t.Pos, OrfNames.GetName(t.Pos), genome)
				}
				unique[genome]++
			}
		}
	}
	for k, v := range unique {
		fmt.Printf("%d has %d unique tags\n", k, v)
	}
}

// Look for any genomes where you share some tags but are 4 muts away on
// others-- this may indicate recombination. Return the genomes you share both
// friends and opposites with.
func FindParadoxes(tags []Tag, which int) []int {
	ret := make([]int, 0)
	friends := make(map[int]bool)
	opposites := make(map[int]bool)

	for _, tag := range tags {
		us := []byte(tag.Mapping[which].Us)
		for k, td := range tag.Mapping {
			nm := utils.NumMuts(us, []byte(td.Us))

			if nm == 0 {
				friends[k] = true
			} else if nm >= 4 {
				opposites[k] = true
			}
		}
	}

	// Now you want the intersection of friends and opposites
	for f, _ := range friends {
		if opposites[f] {
			ret = append(ret, f)
		}
	}

	return ret
}

func ShowAllParadoxes(g *genomes.Genomes, tags []Tag) {
	for i := 0; i < g.NumGenomes(); i++ {
		p := FindParadoxes(tags, i)
		fmt.Printf("%d %d\n", i, len(p))
	}
}

type ParadoxDetail struct {
	pos     int
	other   int
	agrees  bool
	orfName string
}

// paradoxes are the indexes of other genomes with which we have paradoxical
// relationships
func ParadoxDetails(g *genomes.Genomes,
	tags []Tag, which int, paradoxes []int) []ParadoxDetail {

	ret := make([]ParadoxDetail, 0)
	sequenceSimilarities := make(map[int]float64)

	var getSS func(int) float64
	getSS = func(other int) float64 {
		ss, there := sequenceSimilarities[other]
		if !there {
			sequenceSimilarities[other] =
				g.SequenceSimilarity(which, other) * 100
			return getSS(other)
		}
		return ss
	}

	for _, p := range paradoxes {
		for _, tag := range tags {
			us, there := tag.Mapping[which]
			if !there {
				continue
			}

			them, there := tag.Mapping[p]
			if !there {
				continue
			}

			nm := utils.NumMuts([]byte(us.Us), []byte(them.Us))
			name := OrfNames.GetName(tag.Pos)
			ss := getSS(p)

			if nm == 0 {
				fmt.Printf("%d: %d(%s) agrees with %d (%.2f%% ss)\n",
					which, tag.Pos, name, p, ss)
				ret = append(ret, ParadoxDetail{tag.Pos, p, true, name})
			} else if nm >= 4 {
				fmt.Printf("%d: %d(%s) opposes %d (%.2f%% ss)\n",
					which, tag.Pos, name, p, ss)
				ret = append(ret, ParadoxDetail{tag.Pos, p, false, name})
			}
		}
	}
	return ret
}

// Look for differences in the agree/oppose counts organized by OrfName
func SpikeSwap(g *genomes.Genomes, which int, details []ParadoxDetail) {
	for i := 0; i < g.NumGenomes(); i++ {
		var agree, oppose int
		var agreeS, opposeS int
		ss := g.SequenceSimilarity(which, i)

		for _, pd := range details {
			if pd.other != i {
				continue
			}

			if pd.orfName == "S" {
				if pd.agrees {
					agreeS++
				} else {
					opposeS++
				}
			} else {
				if pd.agrees {
					agree++
				} else {
					oppose++
				}
			}
		}

		if agree+oppose+agreeS+opposeS == 0 {
			continue
		}

		fmt.Printf("Comparing with %d(%.2f): %d:%d on S %d:%d off-S\n",
			i, ss*100, agreeS, opposeS, agree, oppose)
	}
}

// How many tags between a and b? Also return their sequence similarity
func CountTags(g *genomes.Genomes, a, b int) (int, float64) {
	g2 := g.Filter(a, b)
	numTags := len(FindPatterns(g2, 6, 4))
	ss := g2.SequenceSimilarity(0, 1)
	return numTags, ss
}

// Return the number of simulations actually run
func Simulate(g *genomes.Genomes, a, b int, count int) int {
	numTags, total := 0, 0
	var numSilent int

	for i := 0; i < count; i++ {
		var gm *genomes.Genomes
		gm, numSilent = MakeSimulatedMutant(g, a, b)
		if gm == nil {
			return 0
		}
		patterns := FindPatterns(gm, 6, 4)
		numTags += len(patterns)
		total++

		/*
			tags := CreateTags(gm, patterns)
			for _, tag := range tags {
				tag.Print()
			}
			gm.SaveMulti("test.fasta")
		*/

		/*
			if i % 10 == 0 {
				fmt.Printf("%d/%d\n", i, count)
			}
		*/
	}

	actual, _ := CountTags(g, a, b)
	fmt.Printf("%d vs %d. Average number of tags: %.2f. "+
		"Actual: %d. Silent muts: %d\n",
		a, b,
		float64(numTags)/float64(total),
		actual,
		numSilent)
	return 1
}

func CreateHighlights(patterns []Pattern) []genomes.Highlight {
	ret := make([]genomes.Highlight, len(patterns))
	for i, p := range patterns {
		ret[i] = genomes.Highlight{p.Pos, p.Pos + len(p.Nts), 'v'}
	}
	return ret
}

func FindElsewhere(g *genomes.Genomes, which int, tag *Tag) int {
	var ret int

    // pos := rand.Intn(g.Length()-8) + 1
    pos := tag.Pos

	pat := g.Nts[which][pos-1 : pos+7]
	var s genomes.Search

	for s.Init(g, which, pat, 0.0); !s.End(); s.Next() {
		ret++
	}

	rc := utils.ReverseComplement(pat)
    if !reflect.DeepEqual(rc, pat) {
        for s.Init(g, which, rc, 0.0); !s.End(); s.Next() {
            ret++
        }
    }

	return ret
}

func SelfRecombination(g *genomes.Genomes, iterations int) {
	n := g.NumGenomes()
	var count, total int

	for i := 0; i < iterations; i++ {
		a, b := rand.Intn(n), rand.Intn(n)
		g2 := g.Filter(a, b)
		patterns := FindPatterns(g2, 6, 4)
		tags := CreateTags(g2, patterns)

		for _, t := range tags {
			found := FindElsewhere(g, 0, &t) > 1
			if found {
				count++
			}
			total++
		}
	}
	fmt.Printf("%d/%d %.2f\n", count, total, float64(count)/float64(total))
}

func main() {
	rand.Seed(1) // FIXME this ain't working
	big := true
	save := false
	print := false

	var g *genomes.Genomes

	if big {
		g = genomes.LoadGenomes("../fasta/SARS2-relatives.fasta",
			"../fasta/WH1.orfs", false)

		/*
			g = genomes.LoadGenomes("../fasta/SARS1-relatives.fasta",
				"../fasta/SARS1.orfs", false)
		*/

		/*
			g = genomes.LoadGenomes("../fasta/more_relatives.fasta",
				"../fasta/WH1.orfs", false)
		*/
	} else {
		g = genomes.LoadGenomes("../fasta/relatives.fasta",
			"../fasta/WH1.orfs", false)
	}
	g.RemoveGaps()

	var tags []Tag

	if save {
		patterns := FindPatterns(g, 6, 4)
		tags = CreateTags(g, patterns)

		SaveTags("tags.gob", tags)
		fmt.Printf("Tags saved")
	} else {
		tags = LoadTags(g, "tags.gob")
	}

	/*
			highlights := []genomes.Highlight{
		        {200, 201, 'x'},
				{100, 140, 'v'},
			}

			g.SaveWithTranslation("test.clu", highlights, 0, 7)
			return
	*/

	if print {
		for _, tag := range tags {
			tag.Print()
		}
	}

	SelfRecombination(g, 1000)

	/*
	   for i := 0; i < g.NumGenomes(); i++ {
	       n, ss := CountTags(g, 7, i)
	       fmt.Printf("7 vs %d (%s): %d tags %.2f%% ss\n", i, g.Names[i], n, ss)
	   }
	*/

	/*
	   g2 := g.Filter(7, 35)
	   patterns := FindPatterns(g2, 6, 4)
	   tags = CreateTags(g2, patterns)
	   highlights := CreateHighlights(patterns)
	   g2.SaveWithTranslation("output.clu", highlights, 0, 1)
	*/

	/*
	   var count int
	   for _, t := range tags {
	       t.Print()
	       if FindElsewhere(g2, 0, &t) > 1 {
	           count++
	       }
	   }
	*/

	/*
		bw := Bandwidth(g, 0, 6, false)
		for i, b := range bw {
			if b >= 5 {
				fmt.Println(i, b)
			}
		}
	*/

	// ShowAllParadoxes(g, tags)
	//FindUnique(tags)

	/*
		paradoxes := FindParadoxes(tags, 0)
		details := ParadoxDetails(g, tags, 0, paradoxes)
		SpikeSwap(g, 0, details)
	*/

	/*
		for i := 0; i < 50; {
			n := g.NumGenomes()
			a, b := rand.Intn(n), rand.Intn(n)
			if a == b {
				continue
			}
			i += Simulate(g, a, b, 500)
		}
	*/

	//g.SaveSelected("WH1-RsYN04.fasta", 0, 54)

	/*
		fmt.Printf("SS inside spike: %.2f%%\n", g.SubSequenceSimilarity(0,
			22, 21562, 25385, true) * 100)

		fmt.Printf("SS outside spike: %.2f%%\n", g.SubSequenceSimilarity(0,
			22, 21562, 25385, false) * 100)
	*/

	/*
		for pos := 28996; pos <= 28999; pos++ {
			t := FindTag(tags, pos)
			t.Compare(0)
		}
	*/
}
