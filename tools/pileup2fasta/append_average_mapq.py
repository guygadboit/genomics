# You run this on the pileup file (which you should use -s to generate)
import sys


def average_mapq(line):
	total, count = 0, 0
	for c in line:
		total += ord(c) - 33
		count += 1
	return float(total)/count


def main():
	with open(sys.argv[1]) as fp:
		for line in fp:
			line = line.strip()
			fields = line.split()
			mapq = fields[6]
			print("{} {:.2f}".format(line, average_mapq(mapq)))


if __name__ == "__main__":
	main()
