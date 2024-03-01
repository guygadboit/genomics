import re
from pdb import set_trace as brk

def main():
	with open("all_pairs_sorted.txt") as fp:
		for line in fp:
			line = line.strip()
			p = float(line.split()[-1])
			m = re.search(r'mps=(\d)', line)
			mps = int(m.group(1))
			print(mps)


if __name__ == "__main__":
	main()
