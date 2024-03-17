from collections import OrderedDict
from copy import *


def parse(fname):
	ret = OrderedDict()
	with open(fname) as fp:
		for line in fp:
			line = line.split(maxsplit=1)
			ret[line[0]] = [float(x) for x in line[1].split()]
		return ret


def combine(a, b):
	ret = a.copy()
	for k, v in b.items():
		if k in ret:
			ret[k].extend(v)
		else:
			ret[k] = v
	return ret


def output(d):
	for k, v in d.items():
		while len(v) < 4:
			v.append(0.0)

		print("{} {}".format(k, " ".join([str(x) for x in v])))


def main():
	sc2 = parse("transitions-SARS2.txt")
	sc1 = parse("transitions-SARS1.txt")
	both = combine(sc2, sc1)
	output(both)


if __name__ == "__main__":
	main()

