#!/usr/bin/bash
seqkit mutate -i $(($1-1)):$2 ../fasta/WH1.fasta > 1.fasta
seqkit subseq -r $(($1-20)):$(($1+32)) 1.fasta > wh1-bit.fasta
