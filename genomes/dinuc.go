package genomes

type DinucProfile struct {
	GC  float64 // Proportion of GC
	CpG float64	// OR of CpG
	TpA float64	// OR of TpA
}

func count(nts []byte, pattern []byte) int {
	m := len(pattern)
	var ret int

searching:
	for i := 0; i < len(nts) - m + 1; i++ {
		for j := 0; j < m; j++ {
			if nts[i+j] != pattern[j] {
				continue searching
			}
		}
		ret++
	}
	return ret
}

func CalcProfile(nts []byte) DinucProfile {
	c := count(nts, []byte{'C'})
	a := count(nts, []byte{'A'})
	t := count(nts, []byte{'T'})
	g := count(nts, []byte{'G'})

	cpgCount := count(nts, []byte("CG"))
	tpaCount := count(nts, []byte("TA"))

	total := float64(c + a + t + g)
	gc := float64(c + g) / total

	// The total number of dinucleotides is 1 less than the total number of
	// nts, which is where that -1 comes from. We're looking at the frequency
	// of CpG over the "expected" frequency.
	totalsq := total * total
	cpg := (float64(cpgCount) / (total - 1.0)) / (float64(c * g) / totalsq)
	tpa := (float64(tpaCount) / (total - 1.0)) / (float64(t * a) / totalsq)

	return DinucProfile{gc, cpg, tpa}
}