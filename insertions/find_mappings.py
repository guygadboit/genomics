from format_raw import load
import os


def main():
	root = os.environ["HOME"] + "/doc/covid/LabLeak/FCS/BieniaszSupplements/"
	data = load("{}/SuppDataSet1.csv".format(root))
	mappings = load("{}/SuppDataSet1Coding.csv".format(root), "read_id")

	for i in [7653, 19792, 20968, 2936, 10475]:
		datum = data[i]
# 		print(mappings[datum["read_id"]])
		print(i, datum["insert_sequence"], len(datum["insert_sequence"]))


if __name__ == "__main__":
	main()
