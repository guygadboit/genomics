from pdb import set_trace as brk

def process(i, sra):
	print("echo Fetching {} {}...".format(i, sra))
	print("fasterq-dump {}".format(sra))
	for j, index in enumerate(("WH1-index",)):
		print("echo aligning with {}...".format(index))
		print("python3 align.py -x {} *.fastq".format(index))
		print("samtools sort -O sam output.sam > sorted.sam")
		print("samtools mpileup sorted.sam > pileup")
		print("pileup2fasta -show pileup | gzip -c > {}-{}.txt.gz".format(
			sra, index))
# 		if j == 0:
# 			print('pileup2fasta -match 23595:CTAATTCACGTA pileup '
# 				'&& echo "FOUND in {}"'.format(sra))
	print("rm *.fastq")
	print("df -h .")
	print("echo processed {}\n".format(i+1))


def process2(i, sra, ref, search_range):
	print("echo Fetching {} {}...".format(i, sra))
	print("fasterq-dump {}".format(sra))
	print('search_reads -v -ref {} -r {}'
	   ' *.fastq && echo "FOUND IN {}"'.format(ref, search_range, sra))
	print("rm *.fastq")
	print("df -h .")
	print("echo processed {}\n".format(i+1))


def main():
	with open("cc_list") as fp:
		for i, line in enumerate(fp):
			fields = line.split()
			name = fields[0]
			process(i, name)
# 			process2(i, name, "WH1-dFCS.fasta", "23591:23630")
#			process2(i, name, "WH1-dRandom.fasta", "23481:23540")


if __name__ == "__main__":
	main()
