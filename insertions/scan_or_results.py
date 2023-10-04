import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plot
from pdb import set_trace as brk
from collections import namedtuple, OrderedDict


SILLY = {
		"GTGTGTGTGTGTGTG",
		"GTGTGTGTG",
		"ACACACA",
		"GTGTGTG",
		"GTGTGTGTGTGTG",
		"AACAAACAAACA",
		"GAGGAGGAGGAG",
		"AGGAGAGAGGAG",
		"GTGTGTGTGTG",
		"AGGAGGAGG",
		"TTCTCTCTCTCT",
		}


Count = namedtuple("Count", "num OR OR2")


def parse_field(f):
	numbers = f.split(",")
	return Count(int(numbers[0]), float(numbers[1]), float(numbers[2]))


def parse():
	ret = {}

	with open("./or-results.txt") as fp:
		for i, line in enumerate(fp):
			line = line.strip()
			fields = line.split()

			if i == 0:
				ret_type = namedtuple("Record", " ".join(fields))
				continue

			pat = fields[1]

			if pat in SILLY:
				continue

			# Patterns like TTTTTTTTTT are also considered silly
			if len(set([c for c in pat])) == 1:
				continue

			rec = ret_type(int(fields[0]), pat,
					*[parse_field(f) for f in fields[2:]])

			ret[pat] = rec

	return ret


def graph(results):
	x, y = OrderedDict(), OrderedDict()
	z = []

	for k, v in results.items():
		interesting = False
		for count in v[2:]:
			# if count.OR >= 10 and count.OR2 >= 10 and count.num > 1:
			if count.OR2 >= 10 and count.num > 1:
				interesting = True
				break

		if interesting:
			for f in v._fields[2:]:
				print("{} {} {} ({})".format(v.id, f,
					getattr(v, f).OR2, len(v.pattern)))

				x[v.id] = True
				y[f] = True
				z.append(getattr(v, f).OR2)


	mp.style.use("ggplot")
	fig = plot.figure(figsize=(10, 10))
	ax = fig.add_subplot(111, projection="3d")

	x_labels = list(x.keys())
	y_labels = list(y.keys())

	_x = np.arange(len(x_labels))
	_y = np.arange(len(y_labels))
	_xx, _yy = np.meshgrid(_x, _y)
	x_points, y_points = _xx.ravel(), _yy.ravel()

	ax.bar3d(x_points, y_points,[0.0] * len(z), 1, 1, z, shade=True)

	xticks = [x + 0.5 for x in range(len(x_labels))]
	ax.set_xticks(xticks, minor=False)
	ax.set_xticklabels(x_labels, rotation=45)

	yticks = [y + 0.5 for y in range(len(y_labels))]
	ax.set_yticks(yticks, minor=False)
	ax.set_yticklabels(y_labels, rotation=90)

	ax.set_title('Title')
# 	plot.show()


def main():
	results = parse()
	graph(results)

if __name__ == "__main__":
	main()
