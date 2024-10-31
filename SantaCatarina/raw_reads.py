import gzip
from pdb import set_trace as brk
import xml.etree.ElementTree as ET


def output_sample(sample):
	s = "\n".join(sample)
	try:
		root = ET.fromstring(s)
	except:
		return
	rr_name, epi_name = None, None
	for Id in root.iter("Id"):
		a = Id.attrib
		if a.get("is_primary") == "1":
			rr_name = Id.text
		elif a.get("db_label") == "Sample name":
			epi_name = Id.text

		if epi_name:
			if epi_name.startswith("EPI_ISL"):
				print(epi_name, rr_name)


def find_sras(fname):
	with gzip.open(fname, "rt") as fp:
		sample = []
		for line in fp:
			line = line.strip()
			if line.find("<BioSample") != -1:
				sample = []

			sample.append(line)

			if line.find("</BioSample") != -1:
				if sample:
					output_sample(sample)


def main():
	find_sras("/fs/f/tmp/biosample_set.xml.gz")


if __name__ == "__main__":
	main()

# I think you're looking for this kind of thing:
# 	Ids>
# 	    <Id db="BioSample" is_primary="1">SAMN14483189</Id>
# 		    <Id db_label="Sample name">EPI_ISL_417917</Id>
# 			    <Id db="SRA">SRS6395995</Id>
