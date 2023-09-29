import sys
import re
import numpy as np
from collections import OrderedDict
import sklearn.decomposition as skd
from pdb import set_trace as brk


def parse_file(fname):
	with open(fname) as fp:
		for line in fp:
			yield line.strip()


class Subsequences:
	def __init__(self, fname):
		m = re.match(r'(\w+)-6.txt', fname)
		self.name = m.group(1)
		self.patterns = {}
		for line in parse_file(fname):
			pat, count = line.split()
			pat = pat[:-1]
			self.patterns[pat] = int(count)

	def sort(self):
		items = list(self.patterns.items())
		items.sort(key=lambda i:i[1])

		self.patterns = OrderedDict()
		for k, v in items:
			self.patterns[k] = v

	def normalize(self):
		total = sum(self.patterns.values())
		for k, v in self.patterns.items():
			self.patterns[k] = float(v) / total

	def display(self, fp=sys.stdout):
		for k, v in self.patterns.items():
			print(k, v, file=fp)

	def align(self, other):
		"""Put self in the same key order as other"""
		new_patterns = OrderedDict()
		for k in other.patterns.keys():
			new_patterns[k] = self.patterns.get(k, 0)
		self.patterns = new_patterns

	def as_array(self):
		ret = np.empty((1, 4096), float)
		for i, v in enumerate(self.patterns.values()):
			ret[0, i] = v

		return ret


def write_reduced_data(subsequences, reduced_data):
	with open("pca.txt", "wt") as fp:
		for ss, row in zip(subsequences, reduced_data):
			print(*row, ss.name, file=fp)
	print("Wrote pca.txt")


def pca(subsequences):
	pca = skd.PCA(2)

	data = subsequences[0].as_array().copy()

	for ss in subsequences[1:]:
		data = np.vstack((data, ss.as_array()))

	pca.fit(data)

	print("Explained variance ratio:", pca.explained_variance_ratio_)
	reduced_data = pca.transform(data)
	write_reduced_data(subsequences, reduced_data)

def main():

	fnames = (
			"Human-6.txt",
			"Bat-6.txt",
			"Pangolin-6.txt",
			"Pig-6.txt",
			"Rabbit-6.txt",
			"RaccoonDog-6.txt",
			"Mouse-6.txt",
			"CC-6.txt",
			"DR-6.txt",
			"HI-6.txt",
			"PA-6.txt",
			"Streptomyces-6.txt",
			"Insertions-6.txt",
			"MaybeBac-6.txt",
			"Listeria-6.txt",
			"Ricksettia-6.txt",
			"Salmonella-6.txt",
			"Legionella-6.txt",
			)

	ss = [Subsequences(f) for f in fnames]
	ss[0].normalize()
	ss[0].sort()

	for s in ss[1:]:
		s.normalize()
		s.align(ss[0])

	pca(ss)

# 	with open("h.txt", "wt") as fp:
# 		h.display(fp)
# 
# 	with open("s.txt", "wt") as fp:
# 		s.display(fp)

if __name__ == "__main__":
	main()






