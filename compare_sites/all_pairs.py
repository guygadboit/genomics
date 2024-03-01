import sys
from pdb import set_trace as brk


def main():
	lines = []
	with open(sys.argv[1]) as fp:
		for line in fp:
			line = line.strip()
			p = float(line.split()[-1])
			lines.append((p, line))
	lines.sort()
	for p, l in lines:
		print(l)


if __name__ == "__main__":
	main()
