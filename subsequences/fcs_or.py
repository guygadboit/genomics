from pdb import set_trace as brk


PATTERN = "CTCCTCGGCGGG"

DATA = (
		("PA", 4, 6264404),
		("CC", 3, 4016947),
		("DR", 4, 2648638),
		("Streptomyces", 12, 9119895),
		("Human", 41, 3298430636),
		("RaccoonDog", 82, 2387097583),
		("Bat", 44, 2059799708),
		("Pangolin", 44, 2059799708),
		)

def expected_frequency(name):
	nts = {}
	with open("{}-1.txt".format(name)) as fp:
		for line in fp:
			line = line.strip()
			fields = line.split()
			nt = fields[0][:-1]
			count = int(fields[1])
			nts[nt] = count

	total = sum(nts.values())
	freq = {k: float(v) / total for k, v in nts.items()}

	ret = 1.0
	for nt in PATTERN:
		ret *= freq[nt]

	return ret


def main():
	for name, a, b in DATA:
		freq = float(a) / b
		e = expected_frequency(name)
		print(name, freq/e)


if __name__ == "__main__":
	main()
