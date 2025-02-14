#!/usr/bin/env python3
import numpy as np
from scipy.stats import fisher_exact
import sklearn.decomposition as skd
import socketserver
import os
import re
import sys
import struct
from pdb import set_trace as brk


SOCK_NAME = "/tmp/call_scipy.sock"


def convert_float(s):
	"""s is the hex representation of an IEEE754 double. Return the actual
	value"""
	v = struct.pack('Q', int(s, 16))
	return struct.unpack('d', v)


def unconvert_float(f):
	"""Return a float as an IEEE754 hex string"""
	v = struct.pack('d', f)
	i = struct.unpack('Q', v)[0]
	return hex(i)[2:]


def parse_array(data):
	ret = None
	rows = data.split(';')
	for row in rows:
		items = np.array([convert_float(x) for x in row.split(',')])
		if ret is None:
			ret = items
		else:
			ret = np.hstack((ret, items))
	return ret.transpose()


def encode_array(arr):
	encoded_rows = []
	for row in arr:
		encoded_row = ",".join([unconvert_float(f) for f in row])
		encoded_rows.append(encoded_row)
	return ";".join(encoded_rows)


class Handler(socketserver.StreamRequestHandler):
	def __init__(self, *args, **kwargs):
		self.pat = re.compile(r'^.*CT\[(\d+) (\d+) (\d+) (\d+)\] (.*)$')
		super(Handler, self).__init__(*args, **kwargs)

	def handle(self):
		print("Client connected...")
		while True:
			line = str(self.rfile.readline(), "utf-8")
			line = line.strip()

			if not line: break

			cmd, line = line.split(maxsplit=1)
			if cmd == "fisher":
				m = self.pat.match(line)

				# We seem to get a blank line when the connection is closed
				if not m: break

				a, b, c, d = [int(x) for x in m.groups()[:4]]
				alternative = m.group(5)
				contingency_table = np.array([[a, b], [c, d]], dtype=float)

				OR, p = fisher_exact(contingency_table, alternative=alternative)
				self.wfile.write(bytes("{} {}\n".format(OR, p), 'ascii'))
			elif cmd == "pca":
				components, line = line.split(maxsplit=1)
				components = int(components)
				data = parse_array(line)
				print("PCA on {} matrix".format(data.shape))

				pca = skd.PCA(components)
				pca.fit(data)

				ev = [unconvert_float(f)
					for f in pca.explained_variance_ratio_]
				transformed_data = pca.transform(data)

				self.wfile.write(bytes("{} {} {}".format(*ev,
						   encode_array(transformed_data)), 'ascii'))


def main():
	try:
		os.unlink(SOCK_NAME)
	except:
		pass
	with socketserver.UnixStreamServer(SOCK_NAME, Handler) as server:
		server.serve_forever()


if __name__ == "__main__":
	main()
