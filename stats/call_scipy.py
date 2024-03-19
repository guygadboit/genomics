#!/usr/bin/env python3
import numpy as np
from scipy.stats import fisher_exact
import socketserver
import os
import re
import sys
from pdb import set_trace as brk


SOCK_NAME = "/tmp/call_scipy.sock"


class Handler(socketserver.StreamRequestHandler):
	def __init__(self, *args, **kwargs):
		self.pat = re.compile(r'^.*CT\[(\d+) (\d+) (\d+) (\d+)\]')
		super(Handler, self).__init__(*args, **kwargs)

	def handle(self):
		print("Client connected...")
		while True:
			line = str(self.rfile.readline())
			m = self.pat.match(line)

			# We seem to get a blank line when the connection is closed
			if not m: break

			a, b, c, d = [int(x) for x in m.groups()]
			contingency_table = np.array([[a, b], [c, d]], dtype=float)
			OR, p = fisher_exact(contingency_table)
			self.wfile.write(bytes("{} {}\n".format(OR, p), 'ascii'))


def main():
	try:
		os.unlink(SOCK_NAME)
	except:
		pass
	with socketserver.UnixStreamServer(SOCK_NAME, Handler) as server:
		server.serve_forever()


if __name__ == "__main__":
	main()
