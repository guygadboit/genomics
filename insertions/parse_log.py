from collections import Counter, namedtuple, defaultdict
import re
from copy import *
from pdb import set_trace as brk


Match = namedtuple("Match", "id pat where E organism")


def parse(fname):
	pat_re = re.compile(r'^([ACTG]+) at (\d+).*')
	record_re = re.compile(r'^(\d+) ([0-9\.]+) (.*)$')
	pat, where = None, None
	ret = []
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			m = pat_re.match(line)
			if m:
				pat, where = m.groups()
				continue
			m = record_re.match(line)
			if m:
				match = Match(int(m.group(1)),
					 pat, where, float(m.group(2)),
					 m.group(3))
# 				if match.E < 1.0:
				ret.append(match)
	return ret


def by_org(records):
	ret = defaultdict(list)
	for rec in records:
		ret[rec.organism].append(rec)
	return ret


def by_E(records):
	by_e = copy(records)
	by_e.sort(key=lambda r: r.E)
	for rec in by_e:
		print(rec)


def count(records):
	counter = Counter()
	for rec in records:
		counter[rec.organism] += 1

	for org, count in counter.most_common(20):
		bo = by_org(records)
		print("{} matches in {}".format(count, org))
		for match in bo[org]:
			print(match)


def main():
	records = parse("2026-08-26.log")
# 	count(records)
	print("Total number of records: {}".format(len(records)))
	by_E(records)


if __name__ == "__main__":
	main()
