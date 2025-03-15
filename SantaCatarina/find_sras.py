import bs4
import sys
import urllib3
from pdb import set_trace as brk


def find_sra(wrapper):
	url = "https://www.ncbi.nlm.nih.gov/sra/?term={}".format(wrapper)
	http = urllib3.PoolManager()
	response = http.request("GET", url)
	soup = bs4.BeautifulSoup(response.data, "lxml")
	for a in soup.find_all("a"):
		href = a.attrs.get("href")
		if "trace.ncbi.nlm.nih.gov" in href:
			return a.getText()

def main():
	# OK sometimes there's more than 1. Then you end up on a page with a couple
	# more links. So you need to think about whether you want to chase that.
	for arg in sys.argv[1:]:
		t = find_sra(arg)
		print(t)
	return

	with open("read_info.txt") as fp:
		for i, line in enumerate(fp):
			line = line.strip()
			fields = line.split()
			sra = find_sra(fields[1]) or ""
			print("<{}> -> <{}>".format(fields[1], sra))
			fields.insert(2, sra)
			print(" ".join(fields))
			if i == 10:
				break

if __name__ == "__main__":
	main()
