import numpy as np
from scipy.stats import fisher_exact
from argparse import ArgumentParser
import codon_optimization as codopt
from collections import namedtuple, Counter, defaultdict
from expectations import *
import re, sys
from pdb import set_trace as brk


def parse_records(fname):
	with open(fname) as fp:
		line = next(iter(fp))
		line = line.strip()
		result_type = namedtuple("Result", line)

		line = next(iter(fp))
		line = line.strip()
		types = line.split()

		for line in fp:
			fields = []
			line = line.strip()
			for t, field in zip(types, line.split()):
				if t == "bool":
					field = field.capitalize()
				else:
					field = "'{}'".format(field)
				fields.append(eval("{}({})".format(t, field)))
			yield result_type(*fields)


def format_record(record):
	d = record._asdict()
	fields = []
	for k, v in d.items():
		if k in ("id", "name", "pattern"):
			fields.append(str(v))
		else:
			fields.append("{}={}".format(k, v))
	fields.append("len={}".format(len(record.pattern)))
	return " ".join(fields)


class Filter:
	def apply(self, record):
		"""Return True for keep"""
		raise NotImplementedError


class Silly(Filter):
	def apply(self, record):
		# Anything with more than 6 of the same nt in a row is probably not
		# interesting
		if re.search(r'(G{6}|A{6}|T{6}|C{6})', record.pattern):
			return False

		# Get rid of simple dinucleotide repeats, like TGTGTGTGTGTG.
		if re.search(r'(..)\1{5,}', record.pattern):
			return False

		return True


class BadPos(Filter):
	def apply(self, record):
		if record.pos >= 29870 or record.pos == 0:
			return False
		return True


class NumHere(Filter):
	def __init__(self, minimum=2):
		self.minimum = minimum

	def apply(self, record):
		return record.num_here >= self.minimum


class StrictNumHere(NumHere):
	def apply(self, record):
		return record.strict_num_here >= self.minimum


class Alone(Filter):
	def __init__(self, min_seqs=2):
		self.min_seqs = min_seqs

	def apply(self, record):
		try:
			return record.seqs >= self.min_seqs
		except AttributeError:
			# This filter only applies to matches, not insertions
			return True


class NumHereOrAlone(Filter):
	def __init__(self, minimum=5):
		self.nh = NumHere(minimum)
		self.alone = Alone()

	def apply(self, record):
		return self.nh.apply(record) or self.alone.apply(record)


class Irrelevant(Filter):
	def apply(self, record):
		try:
			# This is a soil bacterium
			return record.name != "Delftia"
		except AttributeError:
			# This filter only applies to matches, not insertions
			return True


class Length(Filter):
	def __init__(self, minLen, maxLen):
		self.minLen = minLen
		self.maxLen = maxLen

	def apply(self, record):
		n = len(record.pattern)
		return self.minLen <= n <= self.maxLen


class Mul3(Filter):
	def apply(self, record):
		return len(record.pattern) % 3 == 0


class Human(Filter):
	def apply(self, record):
		return not record.in_human


class WH1(Filter):
	def apply(self, record):
		return not record.in_wh1


class HighHomologyHuman(Filter):
	def __init__(self):
		self.ids = set()
		with open("high-human.txt") as fp:
			for line in fp:
				line = line.strip()
				self.ids.add(int(line))

	def apply(self, record):
		ret = record.id not in self.ids
		if not ret:
			print("{} is in the high ids".format(record.id))
		return ret


def filter_records(records, filters):
	ret = []
	for rec in records:
		for filt in filters:
			if not filt.apply(rec):
				break
		else:
			ret.append(rec)
	return ret


def count_insertions(insertions, filters):
	human, wh1, unmatched, total = 0, 0, 0, 0
	for ins in insertions:
		if ins.in_human:
			human += 1
		elif ins.in_wh1:
			wh1 += 1
		elif ins.num_matches == 0:
			unmatched += 1
		total += 1

	print("{} in human, {} in wh1, {} unmatched out of total {}\n".format(
		human, wh1, unmatched, total))

	before = len(insertions)
	insertions = filter_records(insertions, filters)
	n = len(insertions)
	print("Before filters: {}. After filters: {}\n".format(before, n))
	return n


def count_records(records):
	c = Counter()
	hh = 0
	hh_set = set()
	for r in records:
		c[r.id] += 1
		if r.forwards_h + r.backwards_h >= 3:
			hh += 1
			hh_set.add(r.id)
	ids = set([r.id for r in records])
	print("{} matches from {} individual insertions. "
		"{} matches and {} insertions are HH".format(
		len(records), len(ids), hh, len(hh_set)))

	for k, v in c.most_common():
		print("{}: {}".format(k, v))


class Count:
	def __init__(self):
		self.total_matches = 0
		self.total_E = 0.0
		self.total_hE = 0.0

	def add(self, record):
		self.total_matches += 1
		self.total_E += record.E

		if record.hE == 1.0:
			self.total_hE += record.E
		else:
			self.total_hE += record.hE

	def score(self):
		# So this counts really how many more times you appeared than expected
		return self.total_matches / self.total_E


def count_per_organism(records):
	m = defaultdict(Count)
	for rec in records:
		m[rec.name].add(rec)

	print("total matches / E")
	items = list(m.items())
	items.sort(key=lambda i:i[1].score(), reverse=True)
	for k, v in items:
		print(k, v.score())


def count_ends(insertions):
	counts = Counter()
	total = 0
	for ins in insertions:
		counts[ins.pattern[0]] += 1
		counts[ins.pattern[-1]] += 1
		total += 2

	print("Nucleotide frequency at the ends of the inserts")
	for k, v in counts.items():
		print(k, v, v/total)


def distribution(insertions):
	insertions = filter_records(insertions,
						 [StrictNumHere(), Alone(), Mul3()])

	count = 0
	total_nts = 0
	for ins in insertions:
		count += ins.pattern.count("CGGCGG")
		total_nts += len(ins.pattern)

	print("{} out of {} contain CGGCGG ({:.2f} per million nts)\n".format(
		count, len(insertions), (float(count)*1e6)/total_nts))

# 	count_ends(insertions)
# 	output_num_here(insertions)


def cod_human_table():
	cod = CodExp()
	human = HumanExp()
	with open("cod-human-table.xml", "w") as fp:
		print("""
<informaltable>
	<tgroup cols="3">
		<thead>
			<row>
				<entry>Min number of matching nucleotides outside the insertion</entry>
				<entry>Occurrences in cod</entry>
				<entry>Occurrences in human</entry>
			</row>
		</thead>
		<tbody>""", file=fp)

		for n in range(3, 7):
			h = human.get(n)
			c = cod.get(n)
			print("""
			<row>
				<entry>{}</entry>
				<entry>{}/{}</entry>
				<entry>{}/{}</entry>
			</row>""".format(n, *c, *h), file=fp)
		print("""</tbody>
	</tgroup>
</informaltable>""", file=fp)
	print("Wrote cod-human-table.xml")


def unique_records(records, attr):
	seen = set()
	for r in records:
		if getattr(r, attr) in seen:
			continue

		if getattr(r, attr) != '-':
			seen.add(getattr(r, attr))

		yield r


def count_homology(records, names, n):
	passes, total = 0, 0
	for r in records:
		if not names or r.name in names:
			if r.forwards_h + r.backwards_h >= n:
				passes += 1
			total += 1
	return passes, total


def reformat_sf(f):
	s = "{:.2g}".format(f)
	return re.sub(r'e[-0]+(\d+)$', '&#xd7;10<superscript>\\1</superscript>', s)


def homology_frequency(heading, records, names):
	# Use the appropriate ExpShuffle for each thing
	shuffles = [AIsraelExpShuffle(),
			   AViscExpShuffle(),
			   ANaeslExpShuffle(),
			   TreponemaExpShuffle(),
			   PorphyromonasExpShuffle(),
			   TForsythExpShuffle(),
			   AActinomExpShuffle()]
	shuffle_map = {s.name.split()[0]: s for s in shuffles}

	proper_names = {
			"ANaesl": "A.Naeslundii",
			"AIsrael": "A.Israelii",
			"AVisc": "A.Viscosus",
			"Treponema": "T.Denticola",
			"Porphyromonas": "P.Gingivalis",
			"TForsyth": "T.Forsythia",
			"AActinom": "A.Actin.",
			}

	xml_fp = open("comparison.xml", "w")

	with open("homology-ORs.txt", "a") as fp:
		for n in range(3, 7):
			print("n={}".format(n))
			passes, total = count_homology(records, names, n)
			a, b = passes, total

			print("{}: {},".format(n, passes))
			print("<row><entry>{}</entry>".format(n), file=xml_fp)

			results, ps = [], []
			for i, exp in enumerate((CodExp(),
						HumanExp(), shuffle_map.get(heading))):
				# for exp in (CodExp(), HumanExp(), ActinomycesExpMC()):
				if exp is None: continue
				c, d = exp.get(n)

				# Turn it into an odds ratio
				b = total - passes
				d -= c

				contingency_table = np.array([[a, b], [c, d]], dtype=float)
				OR, p = fisher_exact(contingency_table, alternative="greater")

				print("{}: {}/{} passes OR={:.2f} p={:.4g}".format(
						exp.name, passes, total, OR, p))
				results.append(OR)
				ps.append(p)

				if i == 0:
					print("<entry>{}/{}</entry>".format(a, a+b), file=xml_fp)
				print("<entry>OR={:.2f}/p={}</entry>".format(
					OR, reformat_sf(p)), file=xml_fp)
			print("</row>", file=xml_fp)

			all_results = " ".join([str(x) for x in results])
			print("{} {} {}".format(
				proper_names.get(heading, heading), n, all_results), file=fp)

			all_ps = " ".join([str(x) for x in ps])
			print("p: {} {} {}".format(heading, n, all_ps), file=fp)
			print()
		xml_fp.close()


def homology_frequency_table(records):
	names = ("ANaesl",)
	exp = CodExp()
	filters = [
			Silly(),
			BadPos(),
			Irrelevant(),
			Length(12, 24),
			WH1(),
			Mul3(),
			]

	with open("hf_sensitivity_table.xml", "w") as fp:
		for n in [3, 6]:
			print("<simpara>Minimum homology length={}:</simpara>".format(n),
				file=fp)
			print("""<para><informaltable>
			<tgroup cols="5">
				<thead>
					<row>
						<entry>Minimum at this pos \ Minimum seqs</entry>
						<entry>0</entry>
						<entry>1</entry>
						<entry>2</entry>
						<entry>3</entry>
					</row>
				</thead><tbody>""", file=fp)
			for min_nh in range(0, 5):
				print("<row><entry>{}</entry>".format(min_nh), file=fp)
				for min_seqs in range(5):
					extra_filters = [Alone(min_seqs), StrictNumHere(min_nh)]
					f = filter_records(records, filters + extra_filters)
					a, b = count_homology(f, names, n)
					c, d = exp.get(n)

					b -= a
					d -= c

					contingency_table = np.array([[a, b], [c, d]], dtype=float)
					OR, p = fisher_exact(contingency_table,
							alternative="greater")
					print("<entry>OR={:.2f} p={:.2g}</entry>".format(OR, p),
							file=fp)
				print("</row>", file=fp)
			print("</tbody></tgroup></informaltable></para>", file=fp)
	print("Wrote hf_sensitivity_table.xml")


def docbook_record(record, fp):
	species = {"ANaesl": "naesl.",
			"AVisc": "visc.",
			"AIsrael": "israel."}.get(record.name)

	if not species:
		return

	good = False
	if record.forwards_h >= 3:
		good = True
	elif record.backwards_h >= 3:
		good = True

	if not good:
		return

	print("""
	<row>
		   <entry><emphasis>{species}</emphasis></entry>
		   <entry>{record.pattern}</entry>
		   <entry>{record.full_match}</entry>
		   <entry>{record.num_here}</entry>
		   <entry>{record.seqs}</entry>
		   <entry>{record.pos}</entry>
		   <entry>{record.forwards_h}</entry>
		   <entry>{record.backwards_h}</entry>
		   <entry>{record.E}</entry>
		   <entry>{record.hE}</entry>
	</row>""".format(species=species, record=record), file=fp)


def output_num_here(records):
	num_here = {}
	total = 0

	for r in records:
		num_here[r.pos-1] = r.num_here

	with open("num_here.txt", "w") as fp:
		for pos in range(1, 29870):
			n = num_here.get(pos, 0)
			total += n
			print(n, file=fp)

	print("Average num_here: {:.2f}".format(float(total)/29903))
	print("Wrote num_here.txt")

	nh_list = list(num_here.items())
	nh_list.sort(key=lambda x:x[1])

	with open("nh.g", "w") as fp:
		print("var InsertionDistro Distribution = Distribution{", file=fp)
		running_total = 0
		for k, v in nh_list:
			if not (1 < k < 29870):
				continue
			if v <= 40:
				continue
			running_total += v
			print("{{{}, {}}},".format(running_total, k), file=fp)
		print("""}
	""", file=fp)
	print("Wrote nh.g")


def total_homology(record):
	return record.forwards_h + record.backwards_h


def codon_optimization(records):
	passes = {}
	for r in records:
		if codopt.score(r.pattern) >= 8:
			h = total_homology(r)
			try:
				existing, eh = passes[r.id]
				if h >= eh:
					passes[r.id] = (r, h)
			except KeyError:
				passes[r.id] = (r, h)

	for k, v in passes.items():
		r, _ = v
		print(r.pattern, total_homology(r), codopt.score(r.pattern))
# 		print("<row><entry>{}</entry><entry>{}</entry><entry>{}</entry></row>".
# 			format(r.pattern, total_homology(r), codopt.score(r.pattern)))



def main():
	ap = ArgumentParser()
	ap.add_argument("-f", "--optional-filters", default="u")
	ap.add_argument("-s", "--sort-by", default="homology")
	ap.add_argument("-m", "--matches", default="matches.txt")
	ap.add_argument("-d", "--docbook", action="store_true")
	ap.add_argument("-w", "--which", default="all", type=str)
	ap.add_argument("-t", "--tables", action="store_true")
	ap.add_argument("-g", "--gte", default=0, type=int)
	args = ap.parse_args()

	insertions = list(parse_records("insertion-data.txt"))
	distribution(insertions)
	output_num_here(insertions)

	records = list(parse_records(args.matches))

	filters = [
			Silly(),
			BadPos(),
			Irrelevant(),
			Length(12, 200),
			WH1(),
			]

	if "u" in args.optional_filters:
		filters.append(NumHere())
	if "s" in args.optional_filters:
		filters.append(StrictNumHere())
	if "h" in args.optional_filters:
# 		filters.append(Human())
		filters.append(HighHomologyHuman())
	if "3" in args.optional_filters:
		filters.append(Mul3())
	if "a" in args.optional_filters:
		filters.append(Alone())
	if "o" in args.optional_filters:
		filters.append(NumHereOrAlone())

	count = count_insertions(insertions, filters)
	if args.tables:
		homology_frequency_table(records)
		cod_human_table()
		return

	records = filter_records(records, filters)
	count_records(records)

	reverse = False
	if args.sort_by == "homology":
		key = lambda r:r.forwards_h + r.backwards_h
		reverse = True
	elif args.sort_by == 'E':
		key = lambda r:r.E*count

	records.sort(key=key, reverse=reverse)

	if args.gte:
		if args.sort_by != "homology":
			print("You need to sort by homology for this", file=sys.stderr)
			sys.exit(1)

		with open("ids.txt", "w") as fp:
			seen = set()
			for r in records:
				if r.forwards_h + r.backwards_h >= args.gte:
					if r.id not in seen:
						print(r.id, file=fp)
						seen.add(r.id)
				else:
					break
		print("Wrote ids.txt")
		return

	if args.docbook:
		with open("matches.xml", "w") as fp:
			for r in unique_records(records, "full_match"):
				docbook_record(r, fp)
			print("Wrote matches.xml".format(count))
	else:
		for r in unique_records(records, "full_match"):
			print(format_record(r), r.E * count)
	print()

	names = None
	if args.which == "all_act":
		names = set(("AVisc", "ANaesl", "AIsrael"))
	elif args.which == "all":
		names = None
	elif args.which in ("ANaesl", "AVisc", "AIsrael"):
		names = set((args.which,))

	heading = re.sub(r'-matches\.txt', '', args.matches)
	homology_frequency(heading, records, names)
	print()

	try:
		count_per_organism(records)
	except ZeroDivisionError:
		print("Doesn't look like you have E values", file=sys.stderr)


if __name__ == "__main__":
	main()
