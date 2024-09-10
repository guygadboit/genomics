from collections import defaultdict


CODON_TABLE = {
		"TTT": 'F',	# Phenylalanine
		"TTC": 'F',

		"TTA": 'L',	# Leucine
		"TTG": 'L',
		"CTT": 'L',
		"CTC": 'L',
		"CTA": 'L',
		"CTG": 'L',

		"ATT": 'I',	# Isoleucine
		"ATC": 'I',
		"ATA": 'I',

		"ATG": 'M',	# Methionine

		"GTT": 'V',	# Valine
		"GTC": 'V',
		"GTA": 'V',
		"GTG": 'V',

		"TCT": 'S',	# Serine
		"TCC": 'S',
		"TCA": 'S',
		"TCG": 'S',

		"CCT": 'P',	# Proline
		"CCC": 'P',
		"CCA": 'P',
		"CCG": 'P',

		"ACT": 'T',	# Threonine
		"ACC": 'T',
		"ACA": 'T',
		"ACG": 'T',

		"GCT": 'A',	# Alanine
		"GCC": 'A',
		"GCA": 'A',
		"GCG": 'A',

		"TAT": 'Y',	# Tyrosine
		"TAC": 'Y',

		"TAA": '*',	# Stop
		"TAG": '*',

		"CAT": 'H',	# Histidine
		"CAC": 'H',

		"CAA": 'Q',	# Glutadine
		"CAG": 'Q',

		"AAT": 'N',	# Asparagine
		"AAC": 'N',

		"AAA": 'K',	# Lysine
		"AAG": 'K',

		"GAT": 'D',	# Aspartic acid
		"GAC": 'D',

		"GAA": 'E',	# Glutamic acid
		"GAG": 'E',

		"TGT": 'C',	# Cysteine
		"TGC": 'C',

		"TGA": '*',	# Stop
		"TGG": 'W',	# Tryptophan

		"CGT": 'R',	# Arginine
		"CGC": 'R',
		"CGA": 'R',
		"CGG": 'R',

		"AGT": 'S',	# Serine
		"AGC": 'S',

		"AGA": 'R',	# Arginine (again)
		"AGG": 'R',

		"GGT": 'G',	# Glycine
		"GGC": 'G',
		"GGA": 'G',
		"GGG": 'G',
}

def _make_reverse_table():
	ret = defaultdict(list)
	for k, v in CODON_TABLE.items():
		ret[v].append(k)
	return ret

REVERSE_CODON_TABLE = _make_reverse_table()

FREQUENCIES = {
		"TTT": 17.6,
		"TCT": 15.2,
		"TAT": 12.2,
		"TGT": 10.6,
		"TTC": 20.3,
		"TCC": 17.7,
		"TAC": 15.3,
		"TGC": 12.6,
		"TTA":  7.7,
		"TCA": 12.2,
		"TAA":  1.0,
		"TGA":  1.6,
		"TTG": 12.9,
		"TCG":  4.4,
		"TAG":  0.8,
		"TGG": 13.2,
		"CTT": 13.2,
		"CCT": 17.5,
		"CAT": 10.9,
		"CGT":  4.5,
		"CTC": 19.6,
		"CCC": 19.8,
		"CAC": 15.1,
		"CGC": 10.4,
		"CTA":  7.2,
		"CCA": 16.9,
		"CAA": 12.3,
		"CGA":  6.2,
		"CTG": 39.6,
		"CCG":  6.9,
		"CAG": 34.2,
		"CGG": 11.4,
		"ATT": 16.0,
		"ACT": 13.1,
		"AAT": 17.0,
		"AGT": 12.1,
		"ATC": 20.8,
		"ACC": 18.9,
		"AAC": 19.1,
		"AGC": 19.5,
		"ATA":  7.5,
		"ACA": 15.1,
		"AAA": 24.4,
		"AGA": 12.2,
		"ATG": 22.0,
		"ACG":  6.1,
		"AAG": 31.9,
		"AGG": 12.0,
		"GTT": 11.0,
		"GCT": 18.4,
		"GAT": 21.8,
		"GGT": 10.8,
		"GTC": 14.5,
		"GCC": 27.7,
		"GAC": 25.1,
		"GGC": 22.2,
		"GTA":  7.1,
		"GCA": 15.8,
		"GAA": 29.0,
		"GGA": 16.5,
		"GTG": 28.1,
		"GCG":  7.4,
		"GAG": 39.6,
		"GGG": 16.5,
}


def _make_score_table():
	ret = {}
	for aa, codons in REVERSE_CODON_TABLE.items():
		s = []
		for c in codons:
			s.append((FREQUENCIES[c], c))
		s.sort()
		for i, c in enumerate(s):
			ret[c[1]] = i
	return ret

SCORE_TABLE = _make_score_table()


def score(s):
	ret = 0
	for i in range(2, len(s)-1, 3):
		codon = s[i:i+3]
		ret += SCORE_TABLE[codon]
	return ret


def main():
	print(score("CTCCTCGGCGGG"))


if __name__ == "__main__":
	main()
