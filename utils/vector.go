package utils

import (
	"math"
)

func VecAdd(v []float64, other []float64) {
	for i, _ := range v {
		v[i] += other[i]
	}
}

func VecDistance(a, b []float64) float64 {
	var sum float64
	for i, v := range a {
		sum += math.Pow((v - b[i]), 2)
	}
	return math.Sqrt(sum)
}
