import sys
import pysam

"""
module to create output vcf header

"""


def main(bam_header, sample_id, version, contigs):
    vcf_header = []

