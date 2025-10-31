package align

import (
	"bufio"
	"genomics/genomes"
	"os"
	"os/exec"
	"path"
)

// Align them with mafft. We assume the first one has the ORFs you want
func Align(g []*genomes.Genomes,
	tmpDir string) (*genomes.Genomes, error) {
	inputName := path.Join(tmpDir, "unaligned.fasta")
	fd, err := os.Create(inputName)
	if err != nil {
		return nil, err
	}
	defer func() {
		fd.Close()
		os.Remove(inputName)
	}()

	fp := bufio.NewWriter(fd)
	for _, g := range g {
		g.Write(fp, g.Names[0], 0)
	}
	fp.Flush()

	cmd := exec.Command("mafft", "--auto", "--quiet", inputName)
	if err != nil {
		return nil, err
	}

	outputName := path.Join(tmpDir, "aligned.fasta")
	fd2, err := os.Create(outputName)
	if err != nil {
		return nil, err
	}
	defer func() {
		fd2.Close()
		os.Remove(outputName)
	}()
	cmd.Stdout = bufio.NewWriter(fd2)

	err = cmd.Run()
	if err != nil {
		return nil, err
	}
	ret := genomes.LoadGenomes(outputName, "", false)
	ret.Orfs = g[0].Orfs

	return ret, nil
}
