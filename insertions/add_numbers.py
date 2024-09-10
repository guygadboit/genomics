import sys


def main():
	with open(sys.argv[1]) as fp:
		for i, line in enumerate(fp):
			line = line.strip()
			print("{} {}".format(i+1, line))


if __name__ == "__main__":
	 main()
