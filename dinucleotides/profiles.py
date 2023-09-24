from collections import OrderedDict
import string
import os
import subprocess as sp
from pdb import set_trace as brk


GENOME_NAMES = (
		"Insertions",
		"Human",
		"HCoVs",
		"Streptomyces",
		"Haemophilus",
# 		"Bat",
# 		"RaccoonDog",
# 		"Pangolin",
# 		"Rabbit",
# 		"Pig",
# 		"Mouse",
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
		self.gc = self.singles["G"] + self.singles["C"]

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

		# Put them in some consistent order but with CG and TA at the start
		keys = list(sorted(self.doubles.keys()))
		keys.remove("CG")
		keys.remove("TA")
		keys = ["CG", "TA"] + keys

		for k in keys:
			v = self.doubles[k]
			expected = self.singles[k[0]] * self.singles[k[1]]
			actual = v
			ret[k] = actual/expected

		self.profile = ret

	def display(self):
		print(self.name)

		print("G+C: {:.2f}".format(self.gc))

		for k, v in self.profile.items():
			if v <= 0.78:
				sig = " -"
			elif v >= 1.3:
				sig = " +"
			else:
				sig = ""
			print("{}: {:.2f}{}".format(k, v, sig))

		print()

	def plot(self):
		with open("{}.dat".format(self.name), "wt") as fp:
			print("G+C {}".format(self.gc), file=fp)
			for k, v in self.profile.items():
				print("{} {}".format(k, v), file=fp)

		with open("plot.gpi") as fp:
			templ = string.Template(fp.read())

		s = templ.substitute(name=self.name)
		with open("tmp.gpi", "wt") as fp:
			fp.write(s)

		sp.run(["gnuplot", "tmp.gpi"])
		return self.name + ".png"


def main():
	fnames = []
	for name in GENOME_NAMES:
		g = Genome(name)
		g.display()
		fnames.append(g.plot())

	try:
		os.unlink("all.png")
	except FileNotFoundError:
		pass

	sp.run("montage -geometry 640 {} all.png".format(
		" ".join(fnames)), shell=True)
	print("Look at all.png")


if __name__ == "__main__":
	main()
