#!/usr/bin/bash
echo ">query" > query.fasta
echo $2 >> query.fasta

blastn -db=/fs/f/genomes/$1/blast/nucl/nt \
	-max_hsps=1 \
	-query=query.fasta \
	-task=blastn-short
