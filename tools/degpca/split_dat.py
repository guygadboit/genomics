import sys
from argparse import *
from pdb import set_trace as brk


class Redirect:
	def __init__(self, outdir, patterns, fname):
		self.outdir = outdir
		self.patterns = patterns
		self.fname = fname
		self.fp = open("{}{}".format(outdir, fname), "w")

	def _match(self, line):
		for pat in self.patterns:
			if pat in line:
				return True
		return False

	def handle(self, line):
		ret = self._match(line)
		if ret:
			print(line, file=self.fp)
		return ret


class RedirectEverything(Redirect):
	def _match(self, line):
		return True


PATTERNS = (
		"Shanxi",
		"Rspp7",
		(("Wuhan", "BANAL", "RaTG"), "SC2-BANAL-RaTG13"),
		(("Rp140346", "Rs8572", "Rs8561", "Rs8586"), "G4"),
		(("ZC45", "ZXC21"), "Z"),
		)


def main():
	ap = ArgumentParser()
	ap.add_argument("-t", "--title", type=str, default="Unknown")
	ap.add_argument("-o", "--outdir", type=str, default="./")
	ap.add_argument("fname", nargs=1)
	args = ap.parse_args()

	if args.outdir[-1] != '/':
		args.outdir += '/'

	redirects = []
	for item in PATTERNS:
		if isinstance(item, str):
			patterns, name = (item,), item
		else:
			patterns, name = item[0], item[1]
		redirects.append(Redirect(args.outdir, patterns,
				"{}.dat".format(name, ".dat")))

	redirects.append(RedirectEverything(args.outdir, (), "Others.dat"))

	with open(args.fname[0]) as fp:
		for line in fp:
			line = line.strip()
			for r in redirects:
				if r.handle(line):
					break

	with open("{}plot_split.gpi".format(args.outdir), "w") as fp:
		names = ['"{}"'.format(r.fname) for r in redirects]
		print('set title "{}"'.format(args.title), file=fp)
		print("plot {}".format(", ".join(names)), file=fp)

	print("Wrote {}plot_split.gpi".format(args.outdir))


if __name__ == "__main__":
	main()
