from collections import Counter


def load():
	ret = Counter()
	with open("snvs") as fp:
		for line in fp:
			line = line.strip()
			mut = line
			# Reverse these because we think our sample is older
			fr, to = mut[-1], mut[0]
			ret["{}->{}".format(fr, to)] += 1
	return ret


def main():
	spectrum = load()
	for k in (
			"A->C",
			"A->G",
			"A->T",
			"C->A",
			"C->G",
			"C->T",
			"G->A",
			"G->C",
			"G->T",
			"T->A",
			"T->C",
			"T->G",):
		print(k, spectrum[k])


if __name__ == "__main__":
	main()
