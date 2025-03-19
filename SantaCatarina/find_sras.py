import bs4
import sys
import urllib3
from pdb import set_trace as brk


ROOT = "https://www.ncbi.nlm.nih.gov"


def find_sra(soup):
	"""Fish the SRA number out of a page"""
	tr = soup.find("tr", {"class": "sra-run-list-header"})
	if not tr: return None

	table = None
	for parent in tr.parents:
		if parent.name == "table":
			table = parent
			break
	else:
		return None

	for a in table.find_all("a"):
		href = a.attrs.get("href")
		if "trace.ncbi.nlm.nih.gov" in href:
			return a.getText()


def search_multi(soup):
	"""Sometimes there are multiple SRAs and then you get to an intermediate
	page with links"""
	divs = soup.find_all("div", {"class": "rprt"})
	if not divs: raise StopIteration

	for div in divs:
		chilren = div.find_all("div", {"class": "rslt"})
		for child in div.find_all("div", {"class": "rslt"}):
			for link in child.find_all("a"):
				href = link.attrs.get("href")
				yield "{}{}".format(ROOT, href)


def find_sras(http, url):
	print(url)
	response = http.request("GET", url)
	soup = bs4.BeautifulSoup(response.data, "lxml")

	sra = find_sra(soup)
	if sra: yield sra

	if not sra:
		for url in search_multi(soup):
			for sra in find_sras(http, url):
				yield sra


def make_url(wrapper_id):
	return "{}/sra/?term={}".format(ROOT, wrapper_id)


def main():
	# OK sometimes there's more than 1. Then you end up on a page with a couple
	# more links. So you need to think about whether you want to chase that.
	http = urllib3.PoolManager()

# 	for arg in sys.argv[1:]:
# 		for t in find_sras(http, make_url(arg)):
# 			print(t)
# 	return

	with open("read_info2.txt", "wt") as output:
		with open("read_info.txt") as fp:
			for i, line in enumerate(fp):
				line = line.strip()
				fields = line.split()
				sras = list(find_sras(http, make_url(fields[1])))

				if sras:
					sras = ",".join(sras)
				else:
					sras = "NOTFOUND"

				fields.insert(2, sras)
				out_line = " ".join(fields)

				print(out_line)
				print(out_line, file=output)

	print("Wrote read_info2.txt")


if __name__ == "__main__":
	main()
