package genomes

import (
	"bufio"
	"fmt"
	"genomics/utils"
	"io"
	"log"
	"os"
	"strings"
)

/*
Represents a collection of aligned genomes (usually two) with one genome
per row. The orfs "belong" to the first one in the set. We also use this
for a single genome.
*/
type Genomes struct {
	Nts   [][]byte
	Names []string
	Orfs  Orfs
}

func NewGenomes(orfs Orfs, numGenomes int) *Genomes {
	return &Genomes{make([][]byte, numGenomes),
		make([]string, numGenomes), orfs}
}

/*
Load genomes, which might be a fasta file containing a single genome, or
one containing a few of them in an alignment. Be a bit careful when working
with alignments since there may be '-' in there. You might want to call
RemoveGaps before going any futher. If merge load everything into one
genome even if there are multiple > lines (mammal genomes tend to be like
that). orfsName can be empty string if you don't have any ORFs which is
fine if you don't plan on doing any translation.
*/
func LoadGenomes(fname string, orfsName string, merge bool) *Genomes {
	var orfs Orfs

	if orfsName != "" {
		orfs = LoadOrfs(orfsName)
	}

	ret := NewGenomes(orfs, 0)

	fp := utils.NewFileReader(fname)
	defer fp.Close()

	currentRow := make([]byte, 0)
loop:
	// If working with huge genomes uncomment this for faster debugging!
	// for i := 0; i < 1000; i++ {
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)

		if strings.HasPrefix(line, ">") {
			/*
				fields := strings.Fields(line[1:])
				name := fields[0]
			*/
			name := line[1:]
			ret.Names = append(ret.Names, name)
			if !merge && len(currentRow) > 0 {
				ret.Nts = append(ret.Nts, currentRow)
				currentRow = make([]byte, 0)
			}
			continue
		}
		line = strings.ToUpper(line)
		currentRow = append(currentRow, []byte(line)...)
	}
	ret.Nts = append(ret.Nts, currentRow)
	return ret
}

func (g *Genomes) Save(name, fname string, which int) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	fmt.Fprintf(fp, ">%s\n", name)

	nts := g.Nts[which]
	utils.Wrap(fp, nts)
	fp.Flush()
	return nil
}

/*
Make one fasta file with all the nt sequences in it and their names.
*/
func (g *Genomes) SaveMulti(fname string) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)

	for i := 0; i < g.NumGenomes(); i++ {
		fmt.Fprintf(fp, ">%s\n", g.Names[i])
		utils.Wrap(fp, g.Nts[i])
	}
	fp.Flush()
	return nil
}

func (g *Genomes) Clone() *Genomes {
	ret := NewGenomes(g.Orfs, g.NumGenomes())
	for i := 0; i < len(g.Nts); i++ {
		ret.Nts[i] = make([]byte, len(g.Nts[i]))
		copy(ret.Nts[i], g.Nts[i])
	}
	return ret
}

/*
Assuming align is aligned with g, add it to g's own nts array, just doing a
shallow copy
*/
func (g *Genomes) Combine(other *Genomes) {
	for i := 0; i < other.NumGenomes(); i++ {
		g.Nts = append(g.Nts, other.Nts[i])
	}
}

/*
These little functions make it a bit easier not to get confused about
which dimension is which.
*/
func (g *Genomes) NumGenomes() int {
	return len(g.Nts)
}

func (g *Genomes) Length() int {
	return len(g.Nts[0])
}

func (g *Genomes) Slice(which, start, end int) []byte {
	nts := g.Nts[which]

	if start < 0 {
		start = 0
	}
	if end > len(nts) {
		end = len(nts)
	}
	return nts[start:end]
}

/*
Whereever there is a gap in the first genome of an alignment, just remove that
column. You will need to do this for comparisons involving translations to work
as the Orfs don't take into account gaps, and neither do Environment
operations. Returns how many we dropped
*/
func (g *Genomes) RemoveGaps() int {
	var ret int
	nts := g.Nts
	newNts := make([][]byte, g.NumGenomes())

	for i := 0; i < g.NumGenomes(); i++ {
		newNts[i] = make([]byte, 0, g.Length())
	}

	for i := 0; i < g.Length(); i++ {
		if nts[0][i] == '-' {
			ret++
			continue
		}

		for j := 0; j < g.NumGenomes(); j++ {
			newNts[j] = append(newNts[j], nts[j][i])
		}
	}

	g.Nts = newNts
	return ret
}

func (g *Genomes) HaveOrfs() bool {
	return len(g.Orfs) > 0
}
