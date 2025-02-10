package stats

import (
	"fmt"
	"log"
	"math"
	"net"
	"strings"
)

type ContingencyTable struct {
	A, B, C, D int
	OR, P      float64
}

type FisherResult struct {
	OR, p float64
}

type FisherAlternative int

const (
	TWO_SIDED FisherAlternative = iota
	LESS
	GREATER
)

type FisherInput struct {
	ct          *ContingencyTable
	alternative FisherAlternative
}

var fisherInput chan *FisherInput
var fisherOutput chan FisherResult

func parseFloats(buf string, n int) ([]float64, string) {
	components := strings.SplitN(buf, " ", n+1)
	ret := make([]float64, len(components)-1)
	for i := 0; i < n; i++ {
		var x uint64
		_, err := fmt.Sscanf(components[i], "%x", &x)
		if err != nil {
			log.Fatal(err)
		}
		ret[i] = math.Float64frombits(x)
	}
	return ret, components[n]
}

func StatsClient() {
	conn, err := net.Dial("unix", "/tmp/call_scipy.sock")
	if err != nil {
		log.Fatalf("Did you start call_scipy.py?")
	}
	defer conn.Close()

	/*
		Currently this buffer needs to be big enough for everything we expect back
		in one shot
	*/
	buf := make([]byte, 1024*1024)

	alternatives := []string{
		"two-sided",
		"less",
		"greater",
	}

	read := func() {
		_, err := conn.Read(buf)
		if err != nil {
			log.Fatal(err)
		}
	}

	for {
		select {
		case fi := <-fisherInput:
			msg := fmt.Sprintf("fisher %s %s\n",
				fi.ct.String(), alternatives[fi.alternative])
			conn.Write([]byte(msg))
			read()
			var result FisherResult
			fmt.Sscanf(string(buf), "%g %g", &result.OR, &result.p)
			fisherOutput <- result
		case pi := <-pcaInput:
			msg := fmt.Sprintf("pca %d %s\n", pi.components,
				pi.EncodeData())
			fmt.Println(msg)
			conn.Write([]byte(msg))
			read()
			s := string(buf)
			variance, s := parseFloats(s, 2)
			result := PCAResult{variance, nil}
			result.DecodeData(s)

			if len(result.ReducedData) != len(pi.data) {
				// This would happen if your buffer was too small (for example)
				log.Fatalf("Received %d reduced rows from %d rows\n",
					len(result.ReducedData), len(pi.data))
			}
			pcaOutput <- result
		}
	}
}

func (ct *ContingencyTable) Init(a, b, c, d int) {
	ct.A = a
	ct.B = b
	ct.C = c
	ct.D = d
}

// We can work out the OR ourselves
func (c *ContingencyTable) CalcOR() float64 {
	aF, bF, cF, dF := float64(c.A), float64(c.B), float64(c.C), float64(c.D)
	c.OR = (aF / bF) / (cF / dF)
	if math.IsNaN(c.OR) || math.IsInf(c.OR, 1) || math.IsInf(c.OR, -1) {
		c.OR = 0
	}
	return c.OR
}

// But if we want to know the p-value we'll get that from our Python server.
func (c *ContingencyTable) FisherExact(alternative FisherAlternative) (float64,
	float64) {
	fisherInput <- &FisherInput{c, alternative}
	result := <-fisherOutput
	c.OR, c.P = result.OR, result.p
	return c.OR, c.P
}

// A string in the format we pass it to scipy
func (ct *ContingencyTable) String() string {
	return fmt.Sprintf("CT[%d %d %d %d]", ct.A, ct.B, ct.C, ct.D)
}

// A more acceptably readable string
func (ct *ContingencyTable) HumanString() string {
	return fmt.Sprintf("(%d/%d) / (%d/%d)", ct.A, ct.B, ct.C, ct.D)
}

type PCAInput struct {
	components int
	data       [][]float64
}

func (p *PCAInput) EncodeData() string {
	var rows, cols int

	rows = len(p.data)
	if rows > 0 {
		cols = len(p.data[0])
	}

	if rows == 0 || cols == 0 {
		log.Fatal("Invalid PCA Data")
	}

	// Strings for the members of a row
	rowS := make([]string, cols)

	// One string for each row of a matrix
	matrixS := make([]string, rows)

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			rowS[j] = fmt.Sprintf("%x", math.Float64bits(p.data[i][j]))
		}
		matrixS[i] = strings.Join(rowS, ",")
	}
	return strings.Join(matrixS, ";")
}

type PCAResult struct {
	VarianceRatio []float64 // one per component
	ReducedData   [][]float64
}

/*
s is a string representing a matrix in the same format we used in EncodeData: ,
separate items in a row and ; separates rows
*/
func (p *PCAResult) DecodeData(s string) {
	rows := strings.Split(s, ";")
	cols := len(strings.Split(rows[0], ","))

	p.ReducedData = make([][]float64, len(rows))
	for i, row := range rows {
		p.ReducedData[i] = make([]float64, cols)
		for j, item := range strings.Split(row, ",") {
			var v uint64
			fmt.Sscanf(item, "%x", &v)
			p.ReducedData[i][j] = math.Float64frombits(v)
		}
	}
}

var pcaInput chan *PCAInput
var pcaOutput chan PCAResult

func PCA(components int, data [][]float64) PCAResult {
	pcaInput <- &PCAInput{components, data}
	result := <-pcaOutput
	return result
}

func init() {
	fisherInput = make(chan *FisherInput)
	fisherOutput = make(chan FisherResult)

	pcaInput = make(chan *PCAInput)
	pcaOutput = make(chan PCAResult)

	go StatsClient()
}
