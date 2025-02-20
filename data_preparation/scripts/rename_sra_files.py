#! usr/bin/python3
#Kudos to https://github.com/jfq3/Miscellaneous-scripts/blob/master/rename_sra_files.py
# rename_sra_files.py
# John Quensen
# 16 July 2019

# Creates a dictionary from a tab-delimited file of SRA ID's and sample names.
# Uses the dictionary to rename fastq files from SRA sequence IDs to sample names.
# Run the script from the directory containing the fastq sequences and the 
# tab-delimited file of SRA fastq file ID's and sample names.
# Files are assumed to be paired reads with names of the form:
# SRR8648700_1.fastq
# SRR8648700_2.fastq
# Lines in the tab-delimited file should look like (no header line):
# SRR8648702	112-20-1
# SRR8648701	112-20-2
# SRR8648700	112-30-1
# SRR8648699	112-30-2

# Usage is: python3 rename_sra_files.py sra_sample_names_file.tsv

import sys
import os

def main(args):
	name_file = sys.argv[1]
	hash={}
	for line in open(name_file, 'r').readlines():
		file_name, sample_name = line.strip().split()
		if(file_name not in hash.keys()):
			hash[file_name.strip()] = [sample_name.strip()]
		else:
			print("Warning: Duplicate file names..")
			sys.exit()
	for old_file_name in os.listdir():
		if old_file_name.endswith(".fastq"):
			sra = old_file_name.split("_")[0]
			#suf = old_file_name.split("_")[1]
			sam = hash.get(sra)[0]
			new_file_name = "".join(str(sam) + "_" + str(suf))
			os.rename(old_file_name, new_file_name)

if __name__ == "__main__":
        usage = "python3 rename_sra_files.py sra_sample_names_file.tsv"
        if len(sys.argv) != 2 :
            print("Incorrect number of arguments.\nUsage is : ", usage)
            sys.exit()
        main(sys.argv[1:])
