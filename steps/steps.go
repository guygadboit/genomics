package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"genomics/comparison"
	"genomics/genomes"
	"genomics/utils"
	"os"
)

type Window struct {
	Size      int
	Threshold int
}

type WindowData struct {
	Window
	Start   int
	Total   int // Total number of muts in this window
	Missing int // Total number "missing" (indels or Ns)
}

func (w *WindowData) Reset(start int) {
	w.Start = start
	w.Total = 0
}

type WindowDatas struct {
	Datas     []WindowData
	TotalSize int
}

func NewWindowDatas(windows []Window) *WindowDatas {
	datas := make([]WindowData, len(windows))
	total := 0
	for i, w := range windows {
		datas[i] = WindowData{Window: w}
		total += w.Size
	}
	return &WindowDatas{datas, total}
}

func (wd *WindowDatas) Copy() *WindowDatas {
	// TODO I guess you could easily save this reallocation by just copying
	// into an array you already have. But the GC probably works fine anyway.
	datas := make([]WindowData, len(wd.Datas))
	copy(datas, wd.Datas)
	return &WindowDatas{datas, wd.TotalSize}
}

func (wd *WindowDatas) SaveFasta(c *comparison.Comparison, fname string) {
	g := genomes.NewGenomes(nil, 2)
	g.Nts[0] = make([]byte, wd.TotalSize)
	g.Nts[1] = make([]byte, wd.TotalSize)
	start := wd.Datas[0].Start
	for i := 0; i < wd.TotalSize; i++ {
		g.Nts[0][i] = c.Genomes.Nts[c.A][start+i]
		g.Nts[1][i] = c.Genomes.Nts[c.B][start+i]
	}

	makeName := func(index int) string {
		end := wd.TotalSize
		return fmt.Sprintf("%s|%d:%d", c.Genomes.Names[index], start+1, end)
	}

	g.Names[0] = makeName(c.A)
	g.Names[1] = makeName(c.B)

	g.SaveMulti(fname)
	fmt.Printf("Wrote %s\n", fname)
}

func (wd *WindowDatas) Overlaps(other *WindowDatas) bool {
	if other == nil {
		return false
	}
	start := wd.Datas[0].Start
	end := start + wd.TotalSize
	otherStart := other.Datas[0].Start
	otherEnd := otherStart + other.TotalSize
	return min(end, otherEnd) > max(start, otherStart)
}

/*
Look for places where the difference between successive windows exceeds
the thresholds
*/
func MatchWindows(c *comparison.Comparison,
	windows []Window, markers string, csv *csv.Writer, silentOnly bool) bool {
	ret := false
	windowDatas := NewWindowDatas(windows)
	datas := windowDatas.Datas

	length := c.Genomes.Length()
	totalMuts := make([]int, length)    // NS+S counts as we go along
	totalMissing := make([]int, length) // counts of Ns or indels

	i := 0
	c.CumulativeCounts(func(count comparison.Count, cc comparison.Count) {
		nts := c.Genomes.Nts
		totalMuts[i] = count.S
		if !silentOnly {
			totalMuts[i] += count.NS
		}
		if count.Ins > 0 || count.Del > 0 ||
			nts[c.A][i] == 'N' || nts[c.B][i] == 'N' {
			totalMissing[i]++
		}
		i++
	})

	// Fill the windows up starting at pos.
	fillWindowDatas := func(pos int) bool {
		currentWd := 0 // the window we're currently filling
		offset := 0    // how much we've already put into windows

		for j := 0; j < windowDatas.TotalSize; j++ {
			if pos+j == length {
				return false
			}
			wi := j - offset // index into the current window
			if wi == datas[currentWd].Size {
				offset += wi
				wi = 0
				currentWd++
				datas[currentWd].Start = pos + j
			}
			datas[currentWd].Total += totalMuts[pos+j]
			datas[currentWd].Missing += totalMissing[pos+j]
		}
		return true
	}

	// Slide the windows along one position assuming they're already filled at
	// pos-1
	slideWindowDatas := func(pos int) bool {
		if pos+windowDatas.TotalSize > length {
			return false
		}
		for i, _ := range datas {
			datas[i].Total -= totalMuts[pos-1]
			datas[i].Total += totalMuts[pos+datas[i].Size-1]

			datas[i].Missing -= totalMissing[pos-1]
			datas[i].Missing += totalMissing[pos+datas[i].Size-1]

			datas[i].Start += 1
			pos += datas[i].Size
		}
		return true
	}

	var prevInteresting *WindowDatas

	for i := 0; i < length; i++ {
		if i > 0 {
			fillWindowDatas = slideWindowDatas
		}
		if fillWindowDatas(i) {
			var average utils.RollingAverage
			interesting := true
			for i := 1; i < len(datas); i++ {
				difference := datas[i].Total - datas[i-1].Total
				average = average.Add(float64(utils.Abs(difference)))
				threshold := datas[i].Threshold

				// Don't allow missing nts to count
				missing := datas[i-1].Missing + datas[i].Missing
				if difference < 0 {
					difference += missing
					if difference > 0 {
						interesting = false
						break
					}
				} else {
					difference -= missing
					if difference < 0 {
						interesting = false
						break
					}
				}

				if utils.Sign(difference) != utils.Sign(threshold) {
					interesting = false
					break
				}

				if utils.Abs(difference) < utils.Abs(threshold) {
					interesting = false
					break
				}
			}
			if interesting {
				ret = true
				if !windowDatas.Overlaps(prevInteresting) {
					name := fmt.Sprintf("%d-%d", c.A, c.B)
					fname := fmt.Sprintf("%s.txt", name)
					fmt.Printf("%.2f average %s\n", average.Mean(), fname)
					for _, wd := range datas {
						fmt.Println(wd)
					}
					fmt.Println()

					wd := datas[1]
					csv.Write([]string{
						c.Genomes.Names[c.A],
						c.Genomes.Names[c.B],
						utils.Itoa(c.A),
						utils.Itoa(c.B),
						utils.Itoa(wd.Start + 1),
						utils.Itoa(wd.Start + wd.Size),
						utils.Itoa(int(average.Total))})

					c.GraphData(fname)
					c.RunGnuplot(fname,
						markers+JumpMarker(wd.Start, wd.Start+wd.Size),
						true)
					os.Remove(fname)
					windowDatas.SaveFasta(c, fmt.Sprintf("%s.fasta", name))
				}
				prevInteresting = windowDatas.Copy()
			}
		}
	}

	return ret
}

func makeWindows(left, middle, right, threshold int) []Window {
	ret := make([]Window, 3)
	ret[0].Size = left

	ret[1].Size = middle
	ret[1].Threshold = threshold

	ret[2].Size = right
	ret[2].Threshold = -threshold
	return ret
}

func main() {
	var (
		threshold    int
		windowSizesS string
		fasta, orfs  string
		markers      bool
		silentOnly   bool
	)

	flag.StringVar(&fasta, "fasta",
		"../fasta/SARS2-relatives-short-names.fasta", "Alignment")
	flag.StringVar(&orfs, "orfs", "../fasta/WH1.orfs", "ORFs")
	flag.StringVar(&windowSizesS, "windows", "2000,500,2000", "Window sizes")
	flag.IntVar(&threshold, "t", 95, "Threshold")
	flag.BoolVar(&markers, "markers", false, "Add markers (for Spike etc.)")
	flag.BoolVar(&silentOnly, "silent", false, "Consider silent muts only")
	flag.Parse()

	g := genomes.LoadGenomes(fasta, orfs, false)
	g.RemoveGaps()

	windowSizes := utils.ParseInts(windowSizesS, ",")
	windows := makeWindows(windowSizes[0],
		windowSizes[1], windowSizes[2], threshold)

	fd, fp := utils.WriteFile("results.csv")
	defer fd.Close()

	csv := csv.NewWriter(fp)
	csv.Write([]string{"Name A", "Name B",
		"Index A", "Index B", "Start", "End", "Extra muts"})

	var markerString string
	if markers {
		markerString = ORFMarkers(g, "S", "ORF8")
	}

	for i := 0; i < g.NumGenomes(); i++ {
		for j := 0; j < i; j++ {
			c := comparison.Compare(g, i, j)

			if MatchWindows(&c, windows, markerString, csv, silentOnly) {
				fmt.Printf("%s (%d) vs %s (%d)\n", g.Names[i], i, g.Names[j], j)
			}
		}
	}
	fp.Flush()
	fmt.Println("Wrote results.csv")
}
