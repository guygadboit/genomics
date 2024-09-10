package main

import (
	"bufio"
	"fmt"
	"genomics/genomes"
	"genomics/stats"
	"genomics/utils"
	"io"
	"log"
	"os"
	"strings"
)

type Classification int

const (
	ANIMAL   Classification = 1 << 0
	BACTERIA                = 1 << 1
)

type Source struct {
	Classification Classification
	Name           string
	Path           string
	Index          string
	Fasta          string
	NtFreq         map[string]float64
	DinFreq        map[string]float64
	Genome         *genomes.Genomes
}

func (s *Source) LoadGenome() {
	if s.Genome != nil {
		return
	}
	s.Genome = genomes.LoadGenomes(s.Fasta, "", true)
}

func (s *Source) Search(nts []byte, tolerance float64) genomes.Search {
	if tolerance == 0.0 {
		return genomes.NewBidiIndexSearch(s.Index, nts)
	}
	s.LoadGenome()
	return genomes.NewBidiLinearSearch(s.Genome, 0, nts, tolerance)
}

func loadFrequency(source Source, numNts int) map[string]float64 {
	ret := make(map[string]float64)

	f := utils.NewFileReader(fmt.Sprintf(
		"../subsequences/output/%s-%d.txt", source.Name, numNts))
	defer f.Close()

	var total float64
loop:
	for {
		line, err := f.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)
		fields := strings.Fields(line)

		nt := strings.TrimRight(fields[0], ":")
		ret[nt] = float64(utils.Atoi(fields[1]))
		total += ret[nt]
	}

	for k, _ := range ret {
		ret[k] /= total
	}

	return ret
}

func GetSources(class Classification) []Source {
	root := "/fs/f/genomes/"

	animals := []Source{
		{ANIMAL, "Human", "human", root + "human/index",
			root + "human/human.fasta.gz", nil, nil, nil},
		{ANIMAL, "Pangolin", "pangolin", root + "pangolin/index",
			root + "pangolin/pangolin.fasta.gz", nil, nil, nil},
		{ANIMAL, "Bat", "bat", root + "bat/index",
			root + "bat/bat.fasta.gz", nil, nil, nil},
		{ANIMAL, "Raccoon Dog", "raccoon_dog", root + "raccoon_dog/index",
			root + "raccoon_dog/raccoon_dog.fasta.gz", nil, nil, nil},
		{ANIMAL, "Cod", "cod", root + "cod/index",
			root + "cod/cod.fasta.gz", nil, nil, nil},
		{ANIMAL, "Rabbit", "rabbit", root + "rabbit/index",
			root + "rabbit/rabbit.fasta.gz", nil, nil, nil},
		{ANIMAL, "Lizard", "lizard", root + "lizard/index",
			root + "lizard/lizard.fasta.gz", nil, nil, nil},
	}

	bacNames := []string{
		/*
			"Delftia", "Legionella",
			"Salmonella", "Ricksettia", "HI",
			"PA", "Listeria", "Streptomyces",
			"StrepPyogenes", "StrepPneum", "Mycoplasma",
			"Brucella", "OT", "RP",
		*/
		/*
			"Delftia", "Brucella", "StrepPyogenes",
			"StrepPneum", "Mycoplasma", "OT",
			"Porphyromonas", "AActinom", "TForsyth",
			"Treponema", "BactFragilis",
		*/
		"AVisc", "ANaesl", "AIsrael",
		"Treponema", "AActinom", "TForsyth", "Porphyromonas",
	}

	bacteria := make([]Source, len(bacNames))
	for i := 0; i < len(bacNames); i++ {
		name := bacNames[i]
		bacteria[i] = Source{BACTERIA, name,
			fmt.Sprintf("bacteria/%s", name),
			fmt.Sprintf("../fasta/bacteria/%s/index", name),
			fmt.Sprintf("../fasta/bacteria/%s/%s.fasta.gz", name, name),
			nil, nil, nil}
	}

	ret := make([]Source, 0)
	if (class & ANIMAL) != 0 {
		ret = append(ret, animals...)
	}
	if (class & BACTERIA) != 0 {
		ret = append(ret, bacteria...)
	}

	/*
		for i := 0; i < len(ret); i++ {
			ret[i].NtFreq = loadFrequency(ret[i], 1)
			ret[i].DinFreq = loadFrequency(ret[i], 2)
		}
	*/

	return ret
}

/*
What is the expected frequency of pat based on the frequencies of either
the individual nts or the dinucleotides of which it is composed?
*/
func expectedFrequency(pat []byte,
	freq map[string]float64, numNts int) float64 {
	ret := 1.0
	for i := 0; i < len(pat)-numNts; i++ {
		ret *= freq[string(pat[i:i+numNts])]
	}
	return ret
}

func Count(ins *Insertion,
	source *Source, tolerance float64) int {
	var count int
	var search genomes.Search
	for search = source.Search(ins.Nts,
		tolerance); !search.End(); search.Next() {
		count++
	}

	return count
}

func CreateOutput(sources []Source) (*os.File, *bufio.Writer) {
	var outName string

	if len(sources) == 1 {
		outName = fmt.Sprintf("%s-matches.txt", sources[0].Name)
	} else {
		outName = "new-matches.txt"
	}

	f, err := os.Create(outName)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)
	OutputMatchHeader(w)
	fmt.Printf("Writing %s\n", outName)
	return f, w
}

type MatchAction int

const (
	APPEND         MatchAction = 1 << 0
	OUTPUT                     = 1 << 1
	BLAST                      = 1 << 2
	BLAST_HOMOLOGY             = 1 << 3
	RANDOMIZE                  = 1 << 4
)

func CountInGenomes(wh1 *genomes.Genomes,
	id *InsertionData, sources []Source, filters []filterFunc,
	tolerance float64, iterations int, actions MatchAction) {
	bc := stats.BlastDefaultConfig()
	var (
		f *os.File
		w *bufio.Writer
	)

	seen := make(map[string]bool)
	var numProcessed int

	if wh1 == nil {
		wh1 = genomes.LoadGenomes("../fasta/WH1.fasta",
			"../fasta/WH1.orfs", false)
	}

	if actions&OUTPUT != 0 {
		f, w = CreateOutput(sources)
		defer f.Close()
	}

	for it := 0; it < iterations; it++ {
		if actions&RANDOMIZE != 0 {
			id.Randomize(filters)
		}
		filterInsertions(id.Insertions, filters, func(ins *Insertion) {
			// Unique them as we go along
			nts := string(ins.Nts)
			there := seen[nts]
			if there {
				return
			}
			seen[nts] = true

			for i, _ := range sources {
				source := &sources[i]

				for search := source.Search(ins.Nts,
					tolerance); !search.End(); search.Next() {
					pos, _ := search.Get()

					var blastResult stats.BlastResult
					if actions&BLAST != 0 {
						results := stats.Blast(bc, source.Path,
							ins.Nts, 1, 1, stats.NOT_VERBOSE)
						if len(results) != 1 {
							continue
						}
						blastResult = results[0]
					}

					homology := FindHomology(wh1, ins,
						source, pos, search.IsForwards(),
						actions&BLAST_HOMOLOGY != 0)

					dir := BACKWARDS
					if search.IsForwards() {
						dir = FORWARDS
					}
					match := Match{blastResult,
						source.Name, pos, dir, homology}

					// This means there was no homology
					if homology.FullMatch == nil {
						homology.E = blastResult.E
					}

					if actions&OUTPUT != 0 {
						match.Output(w, ins)
					}
					if actions&APPEND != 0 {
						ins.Matches = append(ins.Matches, match)
					}
				}
			}

			numProcessed++
			if numProcessed%10 == 0 {
				fmt.Printf("%d Done %d/%d\n",
					it, numProcessed, len(id.Insertions))
			}

		}, false)
	}

	if w != nil {
		w.Flush()
	}
}

func compl(a byte) byte {
	switch a {
	case 'G':
		return 'C'
	case 'C':
		return 'G'
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	}
	return a
}

// Starting at startA in A and startB in B and going forwards or backwards
// return the maximum number of contiguous matches
func homology(a, b []byte, startA, dirA,
	startB, dirB int, complement bool) int {
	var limitA, limitB, ret int

	/*
		fmt.Printf("A: start=%d dir=%d\n", startA, dirA)
		fmt.Printf("B: start=%d dir=%d\n", startB, dirB)
	*/

	switch dirA {
	case 1:
		limitA = len(a)
		if startA >= limitA {
			return 0
		}
	case -1:
		limitA = -1
		if startA <= limitA {
			return 0
		}
	}
	switch dirB {
	case 1:
		limitB = len(b)
		if startB >= limitB {
			return 0
		}
	case -1:
		limitB = -1
		if startB <= limitB {
			return 0
		}
	}

	var isSame func(a, b byte) bool
	if complement {
		isSame = func(a, b byte) bool {
			return a == compl(b)
		}
	} else {
		isSame = func(a, b byte) bool {
			return a == b
		}
	}

scanning:
	for ai, bi := startA,
		startB; ai != limitA && bi != limitB; ai, bi = ai+dirA, bi+dirB {
		if isSame(a[ai], b[bi]) {
			ret++
		} else {
			break scanning
		}
	}
	return ret
}

func GetFullMatch(source *Source, ins *Insertion,
	pos, forwards, backwards int,
	matchIsForwards bool, doBlast bool) (float64, []byte) {
	bc := stats.BlastDefaultConfig()
	n := len(ins.Nts)
	match := make([]byte, n+backwards+forwards)
	copy(match[:backwards],
		source.Genome.Nts[0][pos-backwards:pos])

	if matchIsForwards {
		copy(match[backwards:backwards+n], ins.Nts)
	} else {
		copy(match[backwards:backwards+n], utils.ReverseComplement(ins.Nts))
	}

	copy(match[backwards+n:],
		source.Genome.Nts[0][pos+n:pos+n+forwards])

	if doBlast {
		results := stats.Blast(bc, source.Path, match, 1, 1, stats.NOT_VERBOSE)
		if len(results) == 1 {
			if !matchIsForwards {
				match = utils.ReverseComplement(match)
			}
			return results[0].E, match
		}
		return 1.0, nil
	} else {
		return 1.0, match
	}
}

func FindHomology(wh1 *genomes.Genomes, ins *Insertion,
	source *Source, sourcePos int, matchForwards bool, doBlast bool) Homology {
	h := Homology{0, 0, 1.0, nil}
	source.LoadGenome()
	n := len(ins.Nts)

	insPos := int(ins.Pos) - 1
	// We seem to have some insertions at 0 and some at 29903. That doesn't
	// make any sense (29903 is the actual length). So let's assume 0 just
	// means 1
	if insPos == -1 {
		insPos = 0
	}

	pos := sourcePos
	a, b := source.Genome.Nts[0], wh1.Nts[0]

	if matchForwards {
		h.Backwards = homology(a, b, pos-1, -1,
			insPos-1, -1, false)
		h.Forwards = homology(a, b, pos+n, 1,
			insPos, 1, false)
	} else {
		h.Backwards = homology(a, b, pos-1, -1,
			insPos, 1, true)
		h.Forwards = homology(a, b, pos+n, 1,
			insPos-1, -1, true)
	}

	if h.Forwards+h.Backwards > 0 {
		h.E, h.FullMatch = GetFullMatch(source, ins,
			pos, h.Forwards, h.Backwards, matchForwards, doBlast)
	}

	return h
}

func countSatellites(insertions []Insertion,
	filters []filterFunc, verbose bool) {
	sources := GetSources(BACTERIA)

	fd, err := os.Create("repeat-results.txt")
	if err != nil {
		log.Fatal("Can't create results file")
	}
	w := bufio.NewWriter(fd)

	defer fd.Close()
	filterInsertions(insertions, filters, func(ins *Insertion) {
		for i := 0; i < len(sources); i++ {
			g := genomes.LoadGenomes(sources[i].Fasta, "", true)

			len, count := findSatellites(g, sources[i].Index,
				ins.Nts, sources[i].Name)
			if count > 0 {
				fmt.Fprintf(w, "%s: %d %s %d %d\n", sources[i].Name,
					ins.Id, string(ins.Nts), len, count)
			}

			rc := utils.ReverseComplement(ins.Nts)
			len, count = findSatellites(g, sources[i].Index,
				rc, sources[i].Name)
			if count > 0 {
				fmt.Fprintf(w, "%s: %d (reversed) %s %d %d\n", sources[i].Name,
					ins.Id, string(rc), len, count)
			}
		}
	}, false)
	w.Flush()
}
