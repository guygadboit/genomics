from scipy.stats import pointbiserialr, kstest
from collections import namedtuple, defaultdict
from argparse import ArgumentParser
from pdb import set_trace as brk


def convert_field(s):
	try:
		return {"true": True, "false": False}[s]
	except KeyError:
		try:
			return int(s)
		except ValueError:
			return s


def parse_results(fname):
	with open(fname) as fp:
		while True:
			line = next(iter(fp))
			line = line.strip()
			if line.startswith('#'):
				continue

			result_type = namedtuple("Result", line)
			break

		for line in fp:
			line = line.strip()
			fields = [convert_field(f) for f in line.split()]
			yield result_type(*fields)


def parse_all_results(fname):
	ret = defaultdict(list)
	for res in parse_results(fname):
		ret[res.name].append(res)
	return ret


def score(results, field):
	x, y = [], []
	for result in results:
		x.append(result.acceptable)
		y.append(getattr(result, field))

	pb = pointbiserialr(x, y)
	return pb.correlation, pb.pvalue


def silent_counts(results):
	for k, v in results.items():
		for field in ("muts_in_sites", "total_sites", "total_singles"):
			print(k, field, score(v, field))


def make_graph_files(results):
	for k, v in results.items():
		with open("{}_true.dat".format(k), "wt") as true_fp:
			with open("{}_false.dat".format(k), "wt") as false_fp:
				for result in v:
					fp = true_fp if result.acceptable else false_fp
					print(result.count, result.max_length, file=fp)


def writeOR(name):
	names = ("{}-ORs".format(name), "{}-passing-ORs".format(name))
	with open(names[0], "w") as all_fp:
		with open(names[1], "w") as passing_fp:
			yield
			while True:
				OR, acceptable = yield
				if OR is None:
					break
				if OR == 0.0:
					continue

				print(OR, file=all_fp)
				if acceptable:
					print(OR, file=passing_fp)

	titles = ["All", "Passing"]
	colours = ["blue", "red"]
	for i, sub_name in enumerate(names):
		title = titles[i]
		colour = colours[i]
		with open("{}.gpi".format(sub_name), "w") as fp:
			print("""set boxwidth 0.05 absolute
set style fill solid 1.0 noborder

bin_width = 0.1
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * (bin_number(x) + 0.8)

set title "Distribution of Hotspot ORs on simulated {name} mutants"
set arrow from 1,0 rto graph 0,1 nohead filled lc "red"

# set arrow from 4.07,0 rto graph 0,0.5 backhead filled lc "magenta"
set arrow from 2.85,0 rto graph 0,0.5 backhead filled lc "black"

plot '{sub_name}' using (rounded($1)):(1) smooth \
		frequency with boxes title "{title}" linecolor "{colour}" """.format(
		**locals()), file=fp)


def rates(results, max_count=None,
		  exact_count=None, require_not_interleaved=False, ors=False):
	for k, v in results.items():
		if ors: w = writeOR(k); next(w)
		good, total = 0, 0
		for result in v:
			acceptable = result.acceptable

			# The "acceptable" field in the results is redundant but it's still
			# worth working it out in the Go program just so we can print it
			# out as we go along. So we might as well check it here.
			assert acceptable == (result.unique and result.max_length < 8000)

			if max_count is not None:
				acceptable = acceptable and result.count <= max_count

			if exact_count is not None:
				acceptable = acceptable and result.count == exact_count

			if require_not_interleaved:
				acceptable = acceptable and not result.interleaved

			if acceptable:
				good += 1

			if ors:
				w.send((result.OR, acceptable))

			total += 1

		if ors:
			try:
				w.send((None, None))
			except StopIteration:
				pass

		print("{}: {}/{} {:.4g}%".format(k, good,
			total, float(good * 100) / total))


def parse_positions(s):
	return [int(x) for x in s.strip('[]').split(',')]


def normalized_positions(result):
	positions = [float(x) for x in parse_positions(result.positions)]
	return [x / float(result.genome_len) for x in positions]


def ks(positions):
	return kstest(positions, "uniform").pvalue


def mean_abs_diff(positions):
	n = len(positions) + 1
	d = 0
	for i, pos in enumerate(positions):
		d += abs(pos - (i+1) / n)
	return d * 1/(n-1)


def check_results(results, test, ref_positions):
	ref = test(ref_positions)
	greater = 0
	total = 0.0

	print("Reference value: {:.3f}".format(ref))
	print("Starting point, mean value, % greater than reference")

	for k, v in results.items():
		greater = 0
		total = 0.0
		for result in v:
			positions = normalized_positions(result)
			val = test(positions)
			total += val
			if val > ref:
				greater += 1
		average = total / len(v)
		print("{} {:.3f} {:.3f}%".format(k, average, (greater * 100) / len(v)))
	print()


def added_removed(results):
	print("Average numbers of sites and averages added and removed")
	for k, v in results.items():
		total_added, total_removed, total_count = 0, 0, 0
		for result in v:
			total_added += result.added
			total_removed += result.removed
			total_count += result.count

		n = len(v)
		added = total_added / n
		removed = total_removed / n
		count = total_count / n
		print("{}: {} muts {:.2f} sites {:.2f} added {:.2f} removed".format(k,
		   v[0].num_muts, count, added, removed))
	print()


def parse_originals():
	ret = {}
	with open("restriction_maps.txt") as fp:
		for line in fp:
			line = line.strip()
			name, sites = line.split()
			name = name[:-1]	# get rid of the colon
			sites = parse_positions(sites)
			ret[name] = sites
	return ret


def extra_sites_in_special(results):
	# The 11 places where sites appear in SC2 and its close relatives. Note
	# this is wrong though, because those positions have moved-- we aren't
	# using aligned genomes!
	special = set([2192, 9750, 10443, 11647,
				17328, 17971, 22921, 22922, 23291, 24101, 24508])
	originals = parse_originals()
	total, count = 0, 0

	for k, result_list in results.items():
		existing = set(originals[k])
		for result in result_list:
			positions = set(parse_positions(result.positions))
			added = positions - existing
			outside = added - special
			total += len(outside)
			count += 1
		print("{}: {:.2f} sites added on average outside the "
			"11 locations".format(k, float(total)/count))


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument("-s", "--silent-counts", action="store_true")
	ap.add_argument("-g", "--graph", action="store_true")
	ap.add_argument("-r", "--rates", action="store_true", default=True)
	ap.add_argument("-m", "--max-count", type=int)
	ap.add_argument("-e", "--exact-count", type=int)
	ap.add_argument("-i", "--require-not-interleaved", action="store_true")
	ap.add_argument("-x", "--extras", action="store_true")
	ap.add_argument("--ors", action="store_true")

	args = ap.parse_args()
	results = parse_all_results(args.fname[0])

	if args.extras:
		extra_sites_in_special(results)

	if args.silent_counts:
		silent_counts(results)

	if args.graph:
		make_graph_files(results)

	if args.max_count or args.exact_count or args.require_not_interleaved:
		rates(results, args.max_count,
			args.exact_count, args.require_not_interleaved, args.ors)
	elif args.rates:
		rates(results, args.ors)

	print()
	added_removed(results)

	return

	WH1 = [0.0733036819048256, 0.32605424204929273,
		0.57947363140822, 0.6009764906531118, 0.8059726448851285]

	print("Kolomogorov Smirnov")
	check_results(results, ks, WH1)

	print("Mean Absolute Difference from perfectly uniform")
	check_results(results, mean_abs_diff, WH1)


if __name__ == "__main__":
	main()
