import gzip
import sys
import re
from pdb import set_trace as brk
import xml.etree.ElementTree as ET


def save_sample(sample):
	with open("debug.xml", "wt") as fp:
		for line in sample:
			print(line, file=fp)


def grep(lines, pattern):
	for line in lines:
		if pattern in line:
			return True
	return False


def search_desc(root):
	# Look for the EPI_ISL bit in the Description
	for desc in root.iter("Description"):
		for element in desc.iter():
			if element.text and "EPI_ISL" in element.text:
				return element.text


def search_attrs(root, attr_name):
	# Look for the EPI_ISL bit in the Attributes
	for atts in root.iter("Attributes"):
		for att in atts.iter("Attribute"):
			if attr_name in att.get("attribute_name").lower():
				if att.text and "EPI_ISL" in att.text:
					return att.text


def clean_epi(epi_name):
	pat = re.compile(r'(EPI_ISL_\d+)')
	m = pat.search(epi_name)
	return m.group(1)


def output_sample(sample):
	s = "\n".join(sample)
	try:
		root = ET.fromstring(s)
	except:
		return

	interested = False
	for on in root.iter("Organism"):
		a = on.attrib
		if a.get("taxonomy_id") == "2697049":
			interested = True
			break

	if not interested:
		return

	rr_name, epi_name, collection_date, location = None, None, None, None

	for Id in root.iter("Id"):
		a = Id.attrib
		if a.get("db") == "SRA":
			rr_name = Id.text
		elif a.get("is_primary") == "1":
			if rr_name is None: rr_name = Id.text
		elif a.get("db_label") == "Sample name":
			if "EPI_ISL" in Id.text:
				epi_name = Id.text

	for At in root.iter("Attribute"):
		a = At.attrib
		if a.get("attribute_name") == "collection_date":
			collection_date = At.text
		elif a.get("attribute_name") == "geo_loc_name":
			location = At.text

	if epi_name is None:
		epi_name = search_desc(root)

	if epi_name is None:
		for s in ("gisaid", "virus strain",
				"sample name", "strain", "isolate",
				"source_name", "passage_history", "note"):
			epi_name = search_attrs(root, s)
			if epi_name: break

	if epi_name is not None:
		epi_name = clean_epi(epi_name)

# 	if grep(sample, "EPI_ISL") and epi_name is None:
# 		save_sample(sample)
# 		brk()

	print(epi_name, rr_name, collection_date, location)


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
	find_sras(sys.argv[1])


if __name__ == "__main__":
	main()

# I think you're looking for this kind of thing:
# 	Ids>
# 	    <Id db="BioSample" is_primary="1">SAMN14483189</Id>
# 		    <Id db_label="Sample name">EPI_ISL_417917</Id>
# 			    <Id db="SRA">SRS6395995</Id>
