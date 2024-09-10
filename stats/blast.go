package stats

import (
	"bufio"
	"crypto/md5"
    "math/rand"
	"fmt"
	"genomics/utils"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"
)

type BlastConfig struct {
	BinDir string // where your blast binaries are
	TmpDir string // where to put query files
	Prefix string // Path prefix for finding genomes
	Suffix string // Path from genome dir to blast dbname
}

func (c *BlastConfig) Init(binDir, tmpDir, prefix, suffix string) {
	c.BinDir = binDir
	c.TmpDir = tmpDir
	c.Prefix = prefix
	c.Suffix = suffix
}

func BlastDefaultConfig() *BlastConfig {
	var ret BlastConfig
	ret.Init("/fs/f/Downloads/ncbi-blast-2.16.0+/bin",
		"/fs/f/tmp", "/fs/f/genomes", "blast/nucl/nt")
	return &ret
}

type BlastResult struct {
	Organism    string
	PctIdentity float64
	Length      int
	Mismatches  int
	E           float64
	Score       float64
}

type BlastResults []BlastResult

func makeName(query []byte) string {
    buf := make([]byte, len(query))

    rand.Read(buf)

    // Overkill-- the random buffer would be fine
    for i := 0; i < len(buf); i++ {
        buf[i] ^= query[i]
    }

	return fmt.Sprintf("%x", md5.Sum(buf))
}

func writeFasta(c *BlastConfig, query []byte) string {
	name := makeName(query)
	fname := name + ".fasta"
	fd, err := os.Create(path.Join(c.TmpDir, fname))
	if err != nil {
		log.Fatal("Can't create temporary file")
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	fmt.Fprintf(fp, ">%s\n", name)
	fmt.Fprintln(fp, string(query))
	fp.Flush()
	return fname
}

type BlastVerbosity int

const (
	VERBOSE BlastVerbosity = iota
	NOT_VERBOSE
)

/*
maxHSP is the maximum number of "high scoring pairs". I think it basically just
means the maximum number of results you want back. If verbose, print the
commands out and don't delete the temporary files.
*/
func Blast(c *BlastConfig, genome string,
	query []byte, maxE float64, maxHSP int,
	verbosity BlastVerbosity) BlastResults {
	verbose := verbosity == VERBOSE
	tmpName := path.Join(c.TmpDir, writeFasta(c, query))
	args := []string{
		fmt.Sprintf("-db=%s", path.Join(c.Prefix, genome, c.Suffix)),
		fmt.Sprintf("-max_hsps=%d", maxHSP),
		fmt.Sprintf("-evalue=%f", maxE),
		fmt.Sprintf("-query=%s", tmpName),
		"-outfmt=6",
	}

	if len(query) < 50 {
		args = append(args, "-task=blastn-short")
	}

	binary := path.Join(c.BinDir, "blastn")

	printCmd := func() {
		fmt.Printf("%s ", binary)
		for _, arg := range args {
			fmt.Printf("%s ", arg)
		}
		fmt.Printf("\n")
	}
	cmd := exec.Command(binary, args...)

	if verbose {
		printCmd()
	}

	stdout, _ := cmd.StdoutPipe()
	stderr, _ := cmd.StderrPipe()
	cmd.Start()

	utils.ReaderLines(stderr, func(line string, err error) bool {
		fmt.Fprintln(os.Stderr, line)
		return true
	})

	ret := make(BlastResults, 0)
	utils.ReaderLines(stdout, func(line string, err error) bool {
		if err != nil {
			log.Fatal(err)
		}
		fields := strings.Fields(line)
		result := BlastResult{
			fields[1],
			utils.Atof(fields[2]),
			utils.Atoi(fields[3]),
			utils.Atoi(fields[4]),
			utils.Atof(fields[10]),
			utils.Atof(fields[11]),
		}
		ret = append(ret, result)
		return true
	})

	err := cmd.Wait()
	if err != nil {
		printCmd()
		log.Fatal("blast returned error")
	} else {
		if !verbose {
			err = os.Remove(tmpName)
			if err != nil {
				log.Fatal(err)
			}
		}
	}

	return ret
}
