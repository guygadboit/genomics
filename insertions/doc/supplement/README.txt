insertions2.txt
===============

These are the actual insertions downloaded from cov-spectrum. Each line starts
with an arbitrary "id" I have given each one that I use in some of the other
files.

*-matches.txt
=============

These are the matches found between the insertions and each bacterium. The
fields are given on the first line, and their Python types (to assist parsing)
on the second. 

The non-obvious fields are explained here:

forwards            True if the insert was found on the + strand
strict_num_here     The number here than had seqs > 1
in_human            Not set (this would be whether the insert matched in human)
forwards_h          Number of bp of adjacent sequence identity before the insert
backwards_h         Number after the insert
score               BLAST score
E                   BLAST E value
hE                  BLAST E value including the adjacent sequence identity

cggcgg.csv
==========

Occurrences of CGGCGG in the inserts without filtering.

cggcgg-filtered.csv
====================

Occurrences of CGGCGG in the inserts with filtering.

fcs-alternatives.txt
====================

Sequence identity and human codon-optimization scores for alternative encodings
of the FCS.
