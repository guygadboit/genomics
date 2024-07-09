import argparse as ap


def parse(fname):
	ret = []
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			ret.append(float(line))
	return ret


def main():
	parser = ap.ArgumentParser()
	parser.add_argument("fname", nargs=1)
	parser.add_argument("-v", "--value", type=float)
	args = parser.parse_args()
	value = args.value

	data = parse(args.fname[0])
	greater, less = 0, 0

	# You should sort these to be quicker but I don't care
	for datum in data:
		if datum < value:
			less += 1
		elif datum > value:
			greater += 1

	print("{} is greater than {} and smaller than {}\n".format(value,
		less, greater))


if __name__ == "__main__":
	main()
