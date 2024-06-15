import sys

def plot(fname):
	with open(fname) as fp:
		with open("temp.dat", "w") as outfp:
			for i, line in enumerate(fp):
				line = line.strip()
				if i == 0:
					title = line
				else:
					print(line, file=outfp)

	with open("plot.gpi", "w") as fp:
		print("""
set title "{}"

set style fill solid
set boxwidth 0.5
set xtics rotate by -90

plot "temp.dat" using 2:xtic(1) with boxes notitle
""".format(title), file=fp)

	print("Now run gnuplot --persist plot.gpi")


def main():
	plot(sys.argv[1])


if __name__ == "__main__":
	main()
