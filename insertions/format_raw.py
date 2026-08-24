import sys
import csv

def main():
	with open(sys.argv[1]) as fp:
		reader = csv.reader(fp)
		headings = next(reader)
		for i, record in enumerate(reader):
			datum = {}
			for h, d in zip(headings, record):
				datum[h] = d

			print("{} ins_{}:{} (1 seqs)".format(i+1,
										datum["ref_coord_before_insertion"],
										datum["insert_sequence"]))

if __name__ == "__main__":
	main()
