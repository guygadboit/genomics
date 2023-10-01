from pdb import set_trace as brk

SILLY = {
		"GTGTGTGTGTGTGTG",
		"GTGTGTGTG",
		"ACACACA",
		"GTGTGTG",
		"GTGTGTGTGTGTG",
		"AACAAACAAACA",
		"GAGGAGGAGGAG",
		"AGGAGAGAGGAG",
		"GTGTGTGTGTG",
		}

def main():
	with open("or-results.txt") as fp:
		for i, line in enumerate(fp):
			line = line.strip()
			fields = line.split()

			if i == 0:
				headings = fields[1:]
				continue

			pat = fields[0]

			if pat in SILLY:
				continue

			# Patterns like TTTTTTTTTT are also considered silly
			if len(set([c for c in pat])) == 1:
				continue

			counts = fields[1:]
			for i, c in enumerate(counts):
				count, OR = [float(x) for x in c.split(',')]
				if OR >= 20.0:
					print("{}: {}".format(headings[i], line))


if __name__ == "__main__":
	main()
