import re
import sys


def read(fname):
	ret = []
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			m = re.match(r'.*p=(.*)$', line)
			if m:
				ret.append(((float(m.group(1))), line))
	ret.sort()
	return ret


def main():
	combined_p = 0
	count = 0
	total_p = 0
	data = read(sys.argv[1])
	for p, line in data:
		print(line)
		total_p += p
		count += 1
	print(total_p, count, total_p / count)


if __name__ == "__main__":
	main()
