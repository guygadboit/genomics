import sys
import re
from collections import defaultdict
from pdb import set_trace as brk


class Sequence:
	def __init__(self, num_muts, acc_no, country, date):
		self.num_muts = num_muts
		self.acc_no = acc_no
		self.country = country
		self.date = date

	def __repr__(self):
		return "{} {} {} {}".format(self.num_muts,
				self.acc_no, self.country, self.date)


def find_sequences(fname):
	ret = defaultdict(list)
	pat = re.compile(r'^(\d+) (EPI_ISL_\d+) (.*) (2020-.*)$')
	with open(fname) as fp:
		for line in fp:
			m = pat.match(line)
			fields = m.groups()
			num, acc, country, date = (int(fields[0]),
				fields[1], fields[2], fields[3])
			ret[num].append(Sequence(num, acc, country, date))
	return ret


def find_reads():
	ret = set()
	with open("/fs/j/genomes/raw_reads/epi_isl.txt") as fp:
		for line in fp:
			line = line.strip()
			ret.add(line)
	return ret


def main():
	reads = find_reads()
	seqs = find_sequences(sys.argv[1])
	for k in sorted(seqs.keys()):
		for seq in seqs[k]:
			if seq.acc_no in reads:
				print(seq)


if __name__ == "__main__":
	main()
