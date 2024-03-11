from scipy.stats import kstest
import sys
from pdb import set_trace as brk


LENGTH = 29903


def analyse_results(fp):
	for line in fp:
		line = line.strip()
		start, numbers = line.split(':')
		data = [float(n) / LENGTH for n in numbers.split()]
		result = kstest(data, "uniform")
		print("{}: {} {:g}".format(start, result.statistic, result.pvalue))


def main():
	if len(sys.argv) > 1:
		with open(sys.argv[1]) as fp:
			analyse_results(fp)
	else:
		analyse_results(sys.stdin)

if __name__ == "__main__":
	main()
