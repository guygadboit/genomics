#!/usr/bin/bash
echo ">query" > query.fasta
echo $1 >> query.fasta

blastn -db=/fs/f/genomes/viruses/SARS2/blast/nucl/nt \
	-max_hsps=1 \
	-query=query.fasta \
	-task=blastn-short

blastn -db=/fs/f/genomes/human/blast/nucl/nt \
	-max_hsps=1 \
	-query=query.fasta \
	-task=blastn-short
