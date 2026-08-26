import sys
import csv
import os
from collections import OrderedDict
from pdb import set_trace as brk


def print_datum(datum):
	print("{} ins_{}:{} (1 seqs)".format(datum["our_id"],
		  datum["ref_coord_before_insertion"], datum["insert_sequence"]))


def load(fname, keyname="our_id"):
	ret = OrderedDict()
	with open(fname) as fp:
		reader = csv.reader(fp)
		headings = next(reader)
		for i, record in enumerate(reader):
			datum = {}
			for h, d in zip(headings, record):
				datum[h] = d
			datum["our_id"] = i+1
			key = datum[keyname]
			ret[key] = datum
		return ret


def main():
	# data = load(sys.argv[1])
	root = os.environ["HOME"] + "/doc/covid/LabLeak/FCS/BieniaszSupplements/"
	data = load("{}/SuppDataSet1.csv".format(root))
# 	mappings = load("{}/SuppDataSet1Coding.csv".format(root), "read_id")

	for datum in data.values():
# 		mapping = mappings[datum["read_id"]]
# 		if mapping["annotation"] == "not_in_cds":
		if datum["position"] == "NA":
			print_datum(datum)


if __name__ == "__main__":
	main()
