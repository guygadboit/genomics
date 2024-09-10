import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plot
from collections import namedtuple, OrderedDict
from pdb import set_trace as brk


Count = namedtuple("Count", "num E")


def parse_field(f):
	numbers = f.split(",")
	return Count(int(numbers[0]), float(numbers[1]))


def is_dinucleotide_repeat(pattern):
	for i in range(2, len(pattern) - 1, 2):
		if pattern[i-2:i] != pattern[i:i+2]:
			return False
	return True


def is_mononucleotide_repeat(pattern):
	return len(set([c for c in pattern])) == 1

species = []

def parse():
	global species
	ret = {}

	with open("./results.txt") as fp:
		for i, line in enumerate(fp):
			line = line.strip()
			fields = line.split()

			if i == 0:
				ret_type = namedtuple("Record", " ".join(fields))
				species = fields[2:]
				continue

			pat = fields[1]

			if (is_mononucleotide_repeat(pat) or
					is_dinucleotide_repeat(pat)):
				continue

			rec = ret_type(int(fields[0]), pat,
					*[parse_field(f) for f in fields[2:]])

			ret[pat] = rec

	return ret


def summary(results):
	x, y = OrderedDict(), OrderedDict()
	z = []

	for k, v in results.items():
		for f in v._fields[2:]:
			count = getattr(v, f)
			print("{} {} {} {} {}".format(v.id, f, count.num, count.E)

			x[v.id] = True
			y[f] = True

def longest_matches(results):
	by_length = []
	for k, v in results.items():
		for count in v[2:]:
			if count.num > 0:
				by_length.append((len(k), v))

	by_length.sort(reverse=True)
	seen = set()

	for k, v in by_length:
		if v.id in seen: continue
		seen.add(v.id)

		values = []
		print("{} {}: ({} nts)".format(v.id, v.pattern,
			len(v.pattern)), end=" ")

		for s in species:
			count = getattr(v, s)
			num = count.num
			if num and count.OR2 > 3:
				values.append("{}: appears {} times (OR: {:.2f} OR2: {:.2f})".
						format(s, num, count.OR, count.OR2))
		print(", ".join(values))


def main():
	results = parse()
	# summary(results)
	longest_matches(results)

if __name__ == "__main__":
	main()
