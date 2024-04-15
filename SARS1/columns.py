import sys


def main():
	for line in sys.stdin:
		pieces = line.split()
		print("{:6}{} {} {} {} {} {} {}".format(pieces[0],
										pieces[1][0],
										pieces[1][1:11],
										pieces[1][11:13],
										pieces[1][13:19],
										pieces[1][19:28],
										pieces[1][28:37],
										pieces[1][37:]))

if __name__ == "__main__":
	main()
