#!/usr/bin/env python3
import sys, re


def main():
	for line in sys.stdin:
		line = line.strip()
		names = [x.strip() for x in line.split(",")]
		print("<authorgroup>")
		for name in names:
			first, last = name.split(maxsplit=1)
			first = re.sub(r'=', ' ', first)
			last = re.sub(r'=', ' ', last)
			print("""<author>
	<personname>
		<firstname>{}</firstname>
		<surname>{}</surname>
	</personname>
 </author>""".format(first, last))
		print("</authorgroup>")


if __name__ == "__main__":
	main()
