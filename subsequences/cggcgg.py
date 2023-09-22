import re
from pdb import set_trace as brk


def reverse_complement(pattern):
	COMPLEMENTS = {
			'G': 'C',
			'A': 'T',
			'C': 'G',
			'T': 'A'
			}

	ret = [COMPLEMENTS[c] for c in pattern]
	return "".join(reversed(ret))


# All the ways of coding for RR
RR = set((
		"CGGCGG",
		"CGTCGT",
		"CGTCGC",
		"CGTCGA",
		"CGTCGG",
		"CGTAGA",
		"CGTAGG",
		"CGCCGT",
		"CGCCGC",
		"CGCCGA",
		"CGCCGG",
		"CGCAGA",
		"CGCAGG",
		"CGACGT",
		"CGACGC",
		"CGACGA",
		"CGACGG",
		"CGAAGA",
		"CGAAGG",
		"CGGCGT",
		"CGGCGC",
		"CGGCGA",
		"CGGAGA",
		"CGGAGG",
		"AGACGT",
		"AGACGC",
		"AGACGA",
		"AGACGG",
		"AGAAGA",
		"AGAAGG",
		"AGGCGT",
		"AGGCGC",
		"AGGCGA",
		"AGGCGG",
		"AGGAGA",
		"AGGAGG",
		))

# And throw the reverse complements in there too
RR |= set([reverse_complement(p) for p in RR])


def parse(fname):
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			pat, count = line.split()
			count = int(count)
			pat = pat[:-1]	# Lose the trailing :
			yield pat, count


def count(fname):
	yes, no = 0, 0
	interested = set(("CGGCGG", "CCGCCG"))

	for pat, count in parse(fname):
		if pat in interested:
			yes += count
		elif pat in RR:
			no += count

	total = yes + no
	return float(yes * 100) / total


def rank(fname):
	print(fname)
	for pat, count in parse(fname):
		if pat in RR:
			print(pat, count)


def main():
	fnames = ["output/" + n for n in (
		"Bat-6.txt",
		"Human-6.txt",
		"Pangolin-6.txt",
		"RaccoonDog-6.txt",
		"Rabbit-6.txt",
		"Pig-6.txt",
		"Mouse-6.txt",
		"Viruses-6.txt",
		"Streptomyces-6.txt",
		"Insertions-6.txt",
		)]

	with open("results.dat", "wt") as fp:
		for fname in fnames:
			# rank(fname)
			m = re.match(r'^output/(\w+)-.*', fname)
			name = m.group(1)
			print(name, count(fname), file=fp)
	print("Now run plot.gpi")


if __name__ == "__main__":
	main()
