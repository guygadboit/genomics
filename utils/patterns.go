package utils

type PatternIterator struct {
	Limit		int		// the limit of each counter (we go < limit)
	Counters	[]int	// the "digits" on the "odometer"
	atEnd		bool
}

func (p *PatternIterator) Init(limit, numCounters int) {
	p.Limit = limit
	p.Counters = make([]int, numCounters)
}

func (p *PatternIterator) Next() {
	for i, _ := range p.Counters {
		p.Counters[i] += 1
		if p.Counters[i] < p.Limit {
			return
		}
		for j := 0; j <= i; j++ {
			p.Counters[i] = 0
		}
	}
	p.atEnd = true
}

func (p *PatternIterator) End() bool {
	return p.atEnd
}

var nts []byte = []byte{'A', 'G', 'T', 'C'}

type NtIterator struct {
	PatternIterator
}

func (ni *NtIterator) Init(length int) {
	ni.PatternIterator.Init(4, length)
}

func (ni *NtIterator) Get() []byte {
	ret := make([]byte, len(ni.Counters))
	for i, v := range ni.Counters {
		ret[i] = nts[v]
	}
	return ret
}
