from collections import namedtuple, defaultdict
import re


Result = namedtuple("Result", "pattern pos name count")


def parse():
	ret = []
	pat = re.compile(r'^([GATC]{6}) at (\d+) in <(.*?)> found (\d+) times')
	with open("results.txt") as fp:
		for line in fp:
			line = line.strip()
			m = pat.match(line)
			nts, pos, name, count = m.groups()
			pos = int(pos)
			count = int(count)
			ret.append(Result(nts, pos, name, count))
	return ret


def by_virus(results):
	by_virus = defaultdict(list)
	for r in results:
		by_virus[r.name].append(r)

	for k, v in by_virus.items():
		total, count = 0, 0
		for r in v:
			total += r.count
			count += 1
			# print(r)
		average = float(total) / count
		special = total > 1 and count == total
		print("Total {} Average count {:.2f} {} {}".format(
				total, average, "special" if special else "", k))


def main():
	print("Results are currently too big to run this!")
	return
	results = parse()
	by_virus(results)


if __name__ == "__main__":
	main()
