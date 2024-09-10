#!/usr/bin/bash
cd /fs/f/genomes/bacteria/$1
export PATH=$PATH:/fs/f/Downloads/ncbi-blast-2.16.0+/bin
zipped_name=$1.fasta.gz
unzipped_name=${zipped_name%.*}
echo "Unzipping..."
gunzip $zipped_name
mkdir -p blast/nucl
cd blast/nucl
echo "Making the db..."
makeblastdb -in ../../$unzipped_name -out nt -dbtype nucl
cd ../..
echo "Zipping again..."
gzip $unzipped_name
