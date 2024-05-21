import itertools
from selenium import webdriver
from selenium.webdriver.common.by import By
import time
from pdb import set_trace as brk


ROOT = ("https://cov-spectrum.org/explore/World/"
		"AllSamples/AllTimes/variants?nucMutations=")


def load():
	ret = []
	with open("silents") as fp:
		for line in fp:
			line = line.strip()
			ret.append(line)
	return ret


def scan(silents):
	driver = webdriver.Chrome()
	for combo in itertools.combinations(silents, 2):
		combo_name = ",".join(combo)
		query = "%2C".join(combo)
		url = ROOT+query
		driver.get(url)
		time.sleep(5)
		try:
			element = driver.find_element(By.CSS_SELECTOR,
					"div#metric-with-tooltip")
			text = element.text
			num = text.split()[0]
			print("{}: {}".format(combo_name, num))
		except:
			print("{}: ERROR".format(combo_name))


def show(combos):
	for c in combos:
		query = "%2C".join(c)
		url = ROOT+query
		print(url)


def main():
	silents = load()
	# scan(silents)

	# The only pairs we found together-- all late and nowhwere near Brazil.
	show((("A5827C", "T14295A"),
	   ("T14295A", "A17799G"),
	   ("A17799G","G21423A"),
	   ("A20481G","G21423A")))


if __name__ == "__main__":
	main()
