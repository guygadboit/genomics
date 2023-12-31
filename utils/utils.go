package utils

import (
	"bufio"
	"compress/gzip"
	"math/rand"
	"log"
	"os"
	"io"
	"fmt"
	"strconv"
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

func RandomNts(length int) []byte {
	nts := [...]byte{'G', 'A', 'T', 'C'}
	ret := make([]byte, length)

	for i := 0; i < length; i++ {
		ret[i] = nts[rand.Intn(len(nts))]
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

func Atoi(s string) int {
	ret, err := strconv.Atoi(s)
	if err != nil {
		log.Fatal("Bad integer")
	}
	return ret
}

func Wrap(w io.Writer, nts []byte) {
	ll := 60

	var i int
	for i = 0; i < len(nts)-ll; i += ll {
		fmt.Fprintf(w, "%s\n", string(nts[i:i+ll]))
	}
	fmt.Fprintf(w, "%s\n", string(nts[i:]))
}
