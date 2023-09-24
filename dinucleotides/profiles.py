from collections import OrderedDict
from pdb import set_trace as brk


GENOME_NAMES = (
		"Viruses",
		"Bat",
		"Human",
		"RaccoonDog",
		"Pangolin",
		"Rabbit",
		"Streptomyces",
		"Pig",
		"Mouse",
		"Insertions",
		"HCoVs",
		)


def lines(fname):
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			yield line


def parse_fname(fname):
	ret = {}
	total = 0
	for line in lines(fname):
		fields = line.split()
		v = int(fields[1])
		ret[fields[0][:-1]] = v
		total += v

	for k, v in ret.items():
		ret[k] = float(v) / total

	return ret


class Genome:
	def __init__(self, name):
		self.name = name

		self.singles = parse_fname("output/{}-1-nts.txt".format(self.name))
		self.doubles = parse_fname("output/{}-2-nts.txt".format(self.name))
		self.merge()
		self.calc_profile()

	def merge(self):
		"""Sources could be inserted in either direction, so merge the counts
		and choose a "canonical" direction"""
		# Note that CG for example, already reverse-complements to CG. So you
		# don't need to merge that one.
		for dest, src in (
				("CC", "GG"),
				("TT", "AA"),
				("TG", "CA"),
				("AG", "CT"),
				("GT", "AC"),
				("GA", "TC"),
				):
			self.doubles[dest] += self.doubles[src]
			del self.doubles[src]

	def calc_profile(self):
		ret = OrderedDict()

		for k, v in sorted(self.doubles.items()):
			expected = self.singles[k[0]] * self.singles[k[1]]
			actual = v
			ret[k] = actual/expected

		self.profile = ret

	def display(self):
		print(self.name)

		for k, v in self.profile.items():
			if v <= 0.78:
				sig = " -"
			elif v >= 1.3:
				sig = " +"
			else:
				sig = ""
			print("{}: {:.6f}{}".format(k, v, sig))

		print()


def main():
	for name in GENOME_NAMES:
		g = Genome(name)
		g.display()


if __name__ == "__main__":
	main()
