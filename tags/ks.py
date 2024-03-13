from scipy.stats import kstest
import sys
from pdb import set_trace as brk


LENGTH = 29903


def normalize(data):
	"""Assume data are sorted low to high"""
	offset = data[0]
	maximum = data[-1] - offset
	if maximum == 0.0:
		return data
	return [(datum - offset) / maximum for datum in data]


def kstest_subdivided(data, num_chunks):
	"""Subdivide data up into chunks and return the average p value of the
	windows"""
	total, count = 0.0, 0
# 	chunk_size = LENGTH // num_chunks
	for i in range(0, LENGTH, chunk_size):
		subdata = data[i:i+chunk_size]
		if subdata:
			try:
				subdata = normalize(data[i:i+chunk_size])
			except:
				continue
			p = kstest(subdata, "uniform").pvalue
			total += p
			count += 1
	return total/count


def analyse_results(fp):
	for line in fp:
		line = line.strip()
		start, numbers = line.split(':')
		data = [float(n) for n in numbers.split()]
		if data:
			stat, p = kstest(normalize(data), "uniform")
			print("{}: {:g} {:g}".format(start, stat, p))


def main():
	if len(sys.argv) > 1:
		with open(sys.argv[1]) as fp:
			analyse_results(fp)
	else:
		analyse_results(sys.stdin)

if __name__ == "__main__":
	main()
