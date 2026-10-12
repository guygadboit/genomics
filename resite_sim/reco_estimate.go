package main

type Locus struct {
	Pos   int
	Count int
	P     float64
}

func (l *Locus) Init(g *genomes.Genomes, pos int, window int) {
	l.Pos = pos
	for i := 1; i < g.NumGenomes(); i++ {
		h := g.Clone()

	}
}

func FindLoci(g *genomes.Genomes) {
	ret := make([]Locus, 0)

	for i := 0; i < g.NumGenomes(); i++ {
	findLocus:
		for j := 0; j < g.Length(); j++ {
			isLocus := -1
			for _, site := range hotspots.RE_SITES {
				for l := 0; l < len(site); l++ {
					if j+l < g.Length() && g.Nts[i][j+l] == site[l] {
						isLocus = j
						break findLocus
					}
				}
			}
			if isLocus != -1 {
				var locus Locus
				locus.Init(g, pos)
				ret = append(ret, locus)
			}
		}
	}
	return ret
}
