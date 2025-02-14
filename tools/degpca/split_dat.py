import sys
from argparse import *
from pdb import set_trace as brk


class Redirect:
	def __init__(self, patterns, fname):
		self.patterns = patterns
		self.fname = fname
		self.writer = self._write(fname); next(self.writer)

	def _match(self, line):
		for pat in self.patterns:
			if pat in line:
				return True
		return False

	def _write(self, fname):
		with open(fname, "w") as fp:
			yield
			while True:
				line = yield
				if line is None:
					break
				print(line, file=fp)

	def handle(self, line):
		ret = self._match(line)
		if ret:
			self.writer.send(line)
		return ret

	def term(self):
		try:
			self.writer.send(None)
		except StopIteration:
			pass


class RedirectEverything(Redirect):
	def _match(self, line):
		return True


PATTERNS = (
		"Shanxi",
		"Rspp7",
		(("Wuhan", "BANAL", "RaTG"), "SC2-BANAL-RaTG13"),
		(("Rp140346", "Rs8572", "Rs8561", "Rs8586"), "G4"),
		"ZC45",
		)


def main():
	ap = ArgumentParser()
	ap.add_argument("-t", "--title", type=str, default="Unknown")
	ap.add_argument("fname", nargs=1)
	args = ap.parse_args()

	redirects = []
	for item in PATTERNS:
		if isinstance(item, str):
			patterns, name = (item,), item
		else:
			patterns, name = item[0], item[1]
		redirects.append(Redirect(patterns, name + ".dat"))

	redirects.append(RedirectEverything((), "Others.dat"))

	with open(args.fname[0]) as fp:
		for line in fp:
			line = line.strip()
			for r in redirects:
				if r.handle(line):
					break

	for r in redirects:
		r.term()

	with open("plot_split.gpi", "w") as fp:
		names = ['"{}"'.format(r.fname) for r in redirects]
		print('set title "{}"'.format(args.title), file=fp)
		print("plot {}".format(", ".join(names)), file=fp)

	print("Wrote plot_split.gpi")


if __name__ == "__main__":
	main()
