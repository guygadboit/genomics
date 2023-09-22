package genomes

import (
	"bufio"
	"compress/gzip"
	"log"
	"os"
	"strings"
)

/*
	Anything iterable and if you do it like this you can easily use it in for
	loops.
*/
type Iter interface {
	Start()                    // Rewind (You probably have your own Init)
	Get() (interface{}, error) // Get the current value
	Next()                     // Advance
	End() bool                 // Are you at the end?
}

func ReverseComplement(nts []byte) []byte {
	ret := make([]byte, len(nts))

	for i := 0; i < len(nts); i++ {

		var nt byte
		switch nts[i] {
		case 'G':
			nt = 'C'
		case 'C':
			nt = 'G'
		case 'A':
			nt = 'T'
		case 'T':
			nt = 'A'
		}

		j := len(nts) - i - 1
		ret[j] = nt
	}
	return ret
}

/*
	Reads files whether they are gzip ones or regular ones
*/
type FileReader struct {
	*bufio.Reader
	fd     *os.File
	gzipFd *gzip.Reader
}

func NewFileReader(fname string) *FileReader {
	var ret FileReader
	ret.Open(fname)
	return &ret
}

func (f *FileReader) Open(fname string) {
	var err error
	f.fd, err = os.Open(fname)
	if err != nil {
		log.Fatalf("Can't open %s", fname)
	}

	if strings.HasSuffix(fname, ".gz") {
		f.gzipFd, err = gzip.NewReader(f.fd)
		if err != nil {
			log.Fatalf("Can't gunzip %s", fname)
		}
		f.Reader = bufio.NewReader(f.gzipFd)
	} else {
		f.Reader = bufio.NewReader(f.fd)
	}
}

func (f *FileReader) Close() {
	if f.gzipFd != nil {
		f.gzipFd.Close()
	}
	f.fd.Close()
}
