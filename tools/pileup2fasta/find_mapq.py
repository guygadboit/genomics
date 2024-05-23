import re
from collections import OrderedDict


def average_mapq(line):
	total, count = 0, 0
	for c in line:
		total += ord(c) - 33
		count += 1
	return float(total)/count


def parse(fp):
	ret = OrderedDict()
	for line in fp:
		line = line.strip()
		m = re.match(r'^[ACTG](\d+)[ACTG]', line)
		pos = int(m.group(1))
		ret[pos] = line
	return ret


def parse_pileup(fp):
	ret = {}
	for line in fp:
		line = line.strip()
		fields = line.split()
		pos = int(fields[1])
		mapq = average_mapq(fields[6])
		ret[pos] = mapq
	return ret


def main():
	with open("muts") as fp:
		muts = parse(fp)
	with open("pileup") as fp:
		qualities = parse_pileup(fp)

	for k, v in muts.items():
		print(v, qualities[k])


if __name__ == "__main__":
	main()
