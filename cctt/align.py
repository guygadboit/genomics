#!/usr/bin/env python3
import subprocess as sp
import os.path
from argparse import *
import sys
from pdb import set_trace as brk


def make_cmd(index, samname, fnames):
	one, two = None, None

	if len(fnames) > 1:
		for fname in fnames:
			base = os.path.splitext(fname)[0]
			if base[-2:] == "_1":
				one = fname
			elif base[-2:] == "_2":
				two = fname

	cmd = ["bowtie2", "--no-unal", "-x", index, "-S", samname]

	if one and two:
		cmd.extend(["-1", one, "-2", two])
	else:
		cmd.append(fnames[0])

	return cmd


def main():
	ap = ArgumentParser()
	ap.add_argument("-x", "--index", type=str)
	ap.add_argument("-s", "--samname", default="output.sam")
	ap.add_argument("fname", type=str, nargs="+")
	args = ap.parse_args()

	cmd = make_cmd(args.index, args.samname, args.fname)
	print(cmd)
	sp.check_call(cmd)


if __name__ == "__main__":
	main()
