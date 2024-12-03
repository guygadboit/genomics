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
	"time"
)

// Let's use this type for OneBased positions and regular ints for zero based.
type OneBasedPos int

// Use this everywhere you need to say if mutations are silent or not
type Silence int
const (
	UNKNOWN Silence = iota
	SILENT
	NON_SILENT
	NOT_IN_ORF
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

func lines(r *bufio.Reader, fun LineFun) {
loop:
	for {
		line, err := r.ReadString('\n')
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

// Call fun for all the lines in a file
func Lines(fname string, fun LineFun) {
	fp := NewFileReader(fname)
	defer fp.Close()
	lines(fp.Reader, fun)
}

// Call fun for all the lines in a Reader
func ReaderLines(reader io.Reader, fun LineFun) {
	r := bufio.NewReader(reader)
	lines(r, fun)
}

func Atoi(s string) int {
	ret, err := strconv.Atoi(s)
	if err != nil {
		log.Fatal("Bad integer")
	}
	return ret
}

func Atof(s string) float64 {
	ret, err := strconv.ParseFloat(s, 64)
	if err != nil {
		log.Fatal("Bad float")
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

// Is it just a mononucleotide repeat?
func IsSilly(pattern []byte) bool {
	for i := 1; i < len(pattern); i++ {
		if pattern[i] != pattern[i-1] {
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

func FromSet[T comparable](s map[T]bool) []T {
	ret := make([]T, 0, len(s))
	for k, _ := range s {
		ret = append(ret, k)
	}
	return ret
}

// a gets b unioned in with it
func Union[T comparable](a map[T]bool, b map[T]bool) {
	for k, _ := range b {
		a[k] = true
	}
}

// Intersection of two sets.
func Intersection[T comparable](a map[T]bool, b map[T]bool) map[T]bool {
	ret := make(map[T]bool)
	for k, _ := range a {
		if b[k] {
			ret[k] = true
		}
	}
	return ret
}

// Return the set a-b
func Difference[T comparable](a map[T]bool, b map[T]bool) map[T]bool {
	ret := make(map[T]bool)
	for k, _ := range a {
		if !b[k] {
			ret[k] = true
		}
	}
	return ret
}

func Shuffle[T any](s []T) {
	rand.Shuffle(len(s), func(i, j int) {
		s[i], s[j] = s[j], s[i]
	})
}

/*
Take a random subset of a slice, but in such a way that if you used the same
random seed you get the same result.
*/
func Sample[T comparable](s []T, count int) []T {
	if count > len(s) {
		return s
	}
	rand.Shuffle(len(s), func(i, j int) {
		s[i], s[j] = s[j], s[i]
	})
	return s[:count]
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

func Shorten(s string, length int) string {
	words := strings.Split(s, " ")
	w := words[0]
	if len(w) > length {
		return w[0:length]
	} else {
		return w
	}
}

/*
Make room in a slice for count new items starting at pos. This is a weird "Go
Idiom". Makes sense when you think about it but too weird to remember so let's
make a function for it.
*/
func Insert[T any](s []T, pos, count int) []T {
	return append(s[:pos], append(make([]T, count), s[pos:]...)...)
}

func Date(year int, month time.Month, day int) time.Time {
	return time.Date(year, month, day, 0, 0, 0, 0, time.UTC)
}

func SplitExt(fname string) (string, string) {
	ext := filepath.Ext(fname)
	base := fname[:len(fname)-len(ext)]
	return base, ext
}

func BaseName(fname string) string {
	var ret, ext string
	for ret = fname; ; {
		ret, ext = SplitExt(ret)
		if len(ext) == 0 {
			return ret
		}
	}
}
