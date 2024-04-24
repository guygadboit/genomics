package utils

type PatternIterator struct {
	limit		int		// the limit of each counter (we go < limit)
	counters	[]int	// the "digits" on the "odometer"
	atEnd		bool
}

func (p *PatternIterator) Init(limit, numCounters int) {
	p.limit = limit
	p.counters = make([]int, numCounters)
}

func (p *PatternIterator) Next() {
	for i, _ := range p.counters {
		p.counters[i] += 1
		if p.counters[i] < p.limit {
			return
		}
		for j := 0; j <= i; j++ {
			p.counters[i] = 0
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
	ret := make([]byte, len(ni.counters))
	for i, v := range ni.counters {
		ret[i] = nts[v]
	}
	return ret
}
