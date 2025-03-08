import bs4
import urllib3
from pdb import set_trace as brk


def find_sra(wrapper):
	url = "https://www.ncbi.nlm.nih.gov/sra/?term={}".format(wrapper)
	http = urllib3.PoolManager()
	response = http.request("GET", url)
	soup = bs4.BeautifulSoup(response.data)
	for a in soup.find_all("a"):
		href = a.attrs.get("href")
		if "trace.ncbi.nlm.nih.gov" in href:
			print(href)
			print(a.attrs.get("href"))
	brk()

def main():
	find_sra("ERS5638945")

if __name__ == "__main__":
	main()
