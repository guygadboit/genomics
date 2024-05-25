# Given the bowtie2 alignment of the reads with WH1 and with P2S create a new
# SAM file containing only the reads that had a higher alignment sore with WH1
# than they did with P2S. Use the non-sorted SAM files for this, so the reads
# are all in the same order
from pdb import set_trace as brk
from collections import namedtuple
from argparse import *
import re


Read = namedtuple("Read", "score line")


def alignment_score(line):
	pat = re.compile(r'^.*AS:i:([-\d]+)')
	m = pat.match(line)

	# No alignment score at all
	if not m: return -9999

	return int(m.group(1))


def load_reads(fname):
	"""Return a dictionary of all the reads"""
	ret = {}

	for line in fname:
		if line.startswith('@'):
			continue
		line = line.strip()
		score = alignment_score(line)
		name = line.split()[0]
		assert name not in ret
		ret[name] = Read(score, line)

	return ret


def find_better(wh1, p2s, output):
	ret = []

	def out(line):
		print(line, file=output)

	sam1_reads = load_reads(wh1)
	sam2_reads = load_reads(p2s)

	out("""@HD\tVN:1.0\tSO:unsorted
@SQ\tSN:Severe\tLN:29903
@PG\tID:me""")

	for k, v in sam1_reads.items():
		other = sam2_reads.get(k)

		if v.score > other.score:
			ret.append((v.score, other.score, v.line))
			out(v.line)

	return ret


def main():
	ap = ArgumentParser()
	ap.add_argument("sam1", nargs=1)
	ap.add_argument("sam2", nargs=1)
	args = ap.parse_args()

	with open(args.sam1[0]) as sam1:
		with open(args.sam2[0]) as sam2:
			with open("better.sam", "wt") as output:
				better = find_better(sam1, sam2, output)
	print("Wrote better.sam")

	# Sort these by how much better they aligned with WH1 than with P2S
	better.sort(key=lambda r: r[0] - r[1], reverse=True)
	for v_score, other_score, line in better:
		print(v_score, other_score, line)


if __name__ == "__main__":
	main()
