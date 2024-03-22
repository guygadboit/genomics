import sys

def main():
	significant = set()
	with open(sys.argv[1]) as fp:
		for line in fp:
			line = line.strip()
			fields = line.split()
			p = float(fields[-1])
			a, b = [int(x) for x in fields[0].split('-')]

			OR = float(fields[-2])
			if OR > 1 and p < 0.05:
				significant |= {a, b}

	if not significant:
		print("None are significantly elevated")
		return

	m = max(significant)
	everything = set([x for x in range(m+1)])
	print(everything - significant)

if __name__ == "__main__":
	main()
