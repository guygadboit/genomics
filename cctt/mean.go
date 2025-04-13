package main

type Mean struct {
	total	float64
	count	int
}

func (m *Mean) Add(value float64) {
	m.total += value
	m.count++
}

func (m *Mean) Get() float64 {
	return m.total / float64(m.count)
}
