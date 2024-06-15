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
	data = read(sys.argv[1])
	for _, line in data:
		print(line)


if __name__ == "__main__":
	main()
