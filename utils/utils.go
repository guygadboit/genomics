package utils

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"path/filepath"
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

type LineFun func(string, error) bool
func Lines(fname string, fun LineFun) {
	fp := NewFileReader(fname)
	defer fp.Close()

loop:
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			line = strings.TrimRight(line, "\n")
			if !fun(line, nil) {
				break
			}
		default:
			fun("", err)
			break
		}
	}
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

// Given a path get the basic filename without the extension out
func BasicName(path string) string {
	_, f := filepath.Split(path)
	ext := filepath.Ext(f)
	return strings.TrimSuffix(f, ext)
}

func IsRegularNt(nt byte) bool {
	switch nt {
	case 'G':
		fallthrough
	case 'A':
		fallthrough
	case 'T':
		fallthrough
	case 'C':
		return true
	default:
		return false
	}
}

func IsRegularPattern(nts []byte) bool {
	for _, nt := range nts {
		if !IsRegularNt(nt) {
			return false
		}
	}
	return true
}

// Given two sequences return the number of differences between them. They are
// assumed to be the same length. Only count "regular" nts, not gaps etc.
func NumMuts(a []byte, b []byte) int {
	var ret int
	for i := 0; i < len(a); i++ {
		if !IsRegularNt(a[i]) || !IsRegularNt(b[i]) {
			continue
		}
		if a[i] != b[i] {
			ret++
		}
	}
	return ret
}

func ToSet[S comparable](s []S) map[S]bool {
	ret := make(map[S]bool)
	for _, item := range s {
		ret[item] = true
	}
	return ret
}

// Parse a , etc. separated list of ints like 0,2,3
func ParseInts(s string, sep string) []int {
	fields := strings.Split(s, sep)
	ret := make([]int, len(fields))
	for i, f := range fields {
		var err error
		ret[i], err = strconv.Atoi(f)
		if err != nil {
			log.Fatalf("Parse error: <%s>\n", s)
		}
	}
	return ret
}
