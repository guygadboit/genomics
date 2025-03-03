from collections import defaultdict


class Sequence:
	def __init__(self, num_muts, acc_no, country):
		self.num_muts = num_muts
		self.acc_no = acc_no
		self.country = country

	def __repr__(self):
		return "{} {} {}".format(self.num_muts, self.acc_no, self.country)


def find_sequences():
	ret = defaultdict(list)
	with open("tt_lots") as fp:
		for line in fp:
			fields = line.strip().split()
			num, acc, country = int(fields[0]), fields[1], fields[2]
			ret[num].append(Sequence(num, acc, country))
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
	seqs = find_sequences()
	for k in sorted(seqs.keys()):
		for seq in seqs[k]:
			if seq.acc_no in reads:
				print(seq)


if __name__ == "__main__":
	main()



