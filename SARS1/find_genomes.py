from collections import namedtuple
import re
from pdb import set_trace as brk


IndexEntry = namedtuple("IndexEntry", "index number")


def parse_index():
	ret = {}
	pat = re.compile(r'^(\d+): ([A-Z_]+\d+).*')
	with open("./index") as fp:
		for line in fp:
			line = line.strip()
			m = pat.match(line)
			index = int(m.group(1))
			number = m.group(2)
			ret[number] = IndexEntry(index, number)
	return ret


def find(index):
	indices = []
	with open("accession_numbers") as fp:
		for line in fp:
			line = line.strip()
			if not line: continue
			number = line.split()[0]
			try:
				indices.append(index[number].index)
			except KeyError:
				print("{} not found".format(number))
	print(",".join([str(x) for x in indices]))


def main():
	index = parse_index()
	find(index)


if __name__ == "__main__":
	main()
