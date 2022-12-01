import sys
import pysam

"""
fetch split reads
sort after split read size
get sequence
compare sequence to input file
"""

def get_region(bamfile, samples, chr, start, end ):
    for read in bamfile.fetch(chr, start, end):
        if read.is_unmapped:
            continue

        if read.is_duplicate:
            continue





def main(bamfile, samples, contigs, contig_length):
    TEs = []

    for chr in contigs: 
        chr_len = contig_length[chr]
        for bin in range(0, chr_len):
            print(bin)

    return TEs
