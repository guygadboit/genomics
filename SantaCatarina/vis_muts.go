package main

import (
	"fmt"
	"genomics/database"
	"genomics/utils"
	"reflect"
	"slices"
)

// Represents a single mutation with a short string code
type MutCode string

// Represents a string of mutations in a unique way
type Code string

type Visualizer struct {
	muts map[string]MutCode
	pi   utils.PatternIterator
}

func initPatternIterator(pi *utils.PatternIterator) {
	pi.Init('z'-'a', 8)
}

func encodePattern(pi utils.PatternIterator) MutCode {
	var ret MutCode = "A"

	for started, i := false, len(pi.Counters)-1; i >= 0; i-- {
		d := pi.Counters[i]

		// Don't bother encoding insignificant zeros-- so any "leading" ones
		// before the actual number starts.
		if d != 0 || started {
			ret += MutCode(byte(d + 'A'))
			started = true
		}
	}
	return ret
}

func (v *Visualizer) Init() {
	v.muts = make(map[string]MutCode)
	initPatternIterator(&v.pi)
}

func (v *Visualizer) Special(m string, c MutCode) {
	v.muts[m] = c
}

func (v *Visualizer) Next() MutCode {
	ret := encodePattern(v.pi)
	v.pi.Next()
	return ret
}

func (v *Visualizer) Encode(muts database.Mutations) Code {
	var ret Code
	codes := make([]MutCode, 0)
	for _, m := range muts {
		ms := m.ToString()
		code, there := v.muts[ms]
		if !there {
			code = v.Next()
			v.muts[ms] = code
		}
		codes = append(codes, code)
	}
	slices.Sort(codes)

	for _, c := range codes {
		ret += Code(fmt.Sprintf("%s-", c))
	}

	return ret
}

func (v *Visualizer) Show() {
	var pi utils.PatternIterator
	initPatternIterator(&pi)

	reverse := make(map[MutCode]string)
	for k, v := range v.muts {
		reverse[v] = k
	}

	for {
		code := encodePattern(pi)
		fmt.Printf("%s: %s\n", code, reverse[code])
		pi.Next()
		if reflect.DeepEqual(pi.Counters, v.pi.Counters) {
			break
		}
	}
}
