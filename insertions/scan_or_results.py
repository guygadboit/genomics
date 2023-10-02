from pdb import set_trace as brk
from collections import namedtuple


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


Count = namedtuple("Count", "num OR OR2")


def parse_field(f):
	numbers = f.split(",")
	return Count(int(numbers[0]), float(numbers[1]), float(numbers[2]))


def parse():
	ret = {}

	with open("./or-results.txt") as fp:
		for i, line in enumerate(fp):
			line = line.strip()
			fields = line.split()

			if i == 0:
				headings = fields[1:]
				ret_type = namedtuple("Record", " ".join(headings))
				continue

			pat = fields[0]

			if pat in SILLY:
				continue

			# Patterns like TTTTTTTTTT are also considered silly
			if len(set([c for c in pat])) == 1:
				continue

			rec = ret_type(*[parse_field(f) for f in fields[1:]])
			ret[pat] = rec

	return ret


def summary(results):
	for k, v in results.items():
		for count in v:
			if count.OR >= 10 and count.OR2 >= 10 and count.num > 1:
				print(k, v, count)

def main():
	results = parse()
	summary(results)

if __name__ == "__main__":
	main()
