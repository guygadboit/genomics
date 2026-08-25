import sys
import csv
from pdb import set_trace as brk


def load(fname, keyname="our_id"):
	ret = {}
	with open(fname) as fp:
		reader = csv.reader(fp)
		headings = next(reader)
		for i, record in enumerate(reader):
			datum = {}
			for h, d in zip(headings, record):
				datum[h] = d

			print("{} ins_{}:{} (1 seqs)".format(i+1,
										datum["ref_coord_before_insertion"],
										datum["insert_sequence"]))
			datum["our_id"] = i+1
			key = datum[keyname]
			ret[key] = datum
		return ret


def main():
	load(sys.argv[1])


if __name__ == "__main__":
	main()
