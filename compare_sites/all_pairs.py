def main():
	lines = []
	with open("all_pairs.txt") as fp:
		for line in fp:
			line = line.strip()
			line.split()
			p = float(line.split()[-1])
			lines.append((p, line))
	lines.sort()
	for p, l in lines:
		print(l)


if __name__ == "__main__":
	main()
