package stats

import (
	"fmt"
	"log"
	"math"
	"net"
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

func FisherClient() {
	conn, err := net.Dial("unix", "/tmp/call_scipy.sock")
	if err != nil {
		log.Fatalf("Did you start call_scipy.py?")
	}
	defer conn.Close()

	buf := make([]byte, 256)
	var result FisherResult

	alternatives := []string{
		"two-sided",
		"less",
		"greater",
	}

	for {
		fi := <-fisherInput
		msg := fmt.Sprintf("%s %s\n",
			fi.ct.String(), alternatives[fi.alternative])
		conn.Write([]byte(msg))
		_, err := conn.Read(buf)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Sscanf(string(buf), "%g %g", &result.OR, &result.p)
		fisherOutput <- result
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

func (ct *ContingencyTable) String() string {
	return fmt.Sprintf("CT[%d %d %d %d]", ct.A, ct.B, ct.C, ct.D)
}

func init() {
	fisherInput = make(chan *FisherInput)
	fisherOutput = make(chan FisherResult)
	go FisherClient()
}
