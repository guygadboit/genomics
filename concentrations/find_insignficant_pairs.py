import sys
from collections import namedtuple, defaultdict
from pdb import set_trace as brk


Record = namedtuple("Record", "a b OR p line")


def parse(fname):
	ret = []
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			fields = line.split()
			a, b = [int(x) for x in fields[0].split('-')]
			p = float(fields[-1])
			OR = float(fields[-2])
			ret.append(Record(a, b, OR, p, line))
	return ret


def find_significant(records):
	significant = set()
	for rec in records:
		if rec.OR > 1 and rec.p < 0.05:
			significant |= {rec.a, rec.b}

	if not significant:
		print("None are significantly elevated")
		return

	m = max(significant)
	everything = set([x for x in range(m+1)])
	print("Insignificant pairs")
	insignificant = everything - significant
	print(insignificant)

	print("Significant pairs")
	print(significant)
	print()


def sort(records):
	"""Show the most significant first"""
	counts = defaultdict(int)
	for rec in records:
		if rec.OR > 1 and rec.p < 0.05:
			counts[rec.a] += 1
			counts[rec.b] += 1

	data = list(counts.items())
	data.sort(key=lambda x: x[1], reverse=True)

	print("Sorted by number of significant pairs involved in:")
	for datum in data:
		print("{}: {}".format(*datum))


def main():
	records = parse(sys.argv[1])
	find_significant(records)
	sort(records)


if __name__ == "__main__":
	main()
