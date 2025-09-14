package comparison

type DataSeries []float64

// Treat the first sample's derivative as 0. This should work quite well if
// you're differentiating cumulative counts.
func (d DataSeries) Differentiate() DataSeries {
	ret := make(DataSeries, len(d))
	for i := 1; i < len(d); i++ {
		ret[i] = d[i] - d[i-1]
	}
	return ret
}
