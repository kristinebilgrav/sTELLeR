import sys
import argparse
import pathlib
import os
import pysam

"""
import py script modules
define main
parsers
setups/wrapper
"""


import tellr.vcf_header as tellr_vcf_header 
import tellr.variants as tellr_call_vars

def main():
    version = "0.0.0"
    parser= argparse.ArgumentParser(
        prog = 'tellr', 
        description='calls non-reference transposable elements given in xx file in long-read pacbio or ont bam files', 
        epilog='some help text'
    )

    parser.add_argument('-ref', '--ref', help='reference genome')
    parser.add_argument('-TE_ref', '--TE_ref', help='fasta sequence of elements to be detected', required=True, type=pathlib.Path)
    parser.add_argument('-bam', '--bam', help='bam file', required= True, type=pathlib.Path)

    args = parser.parse_args()

    if not os.path.isfile(args.bam) :
        print('Error, cannot find bam file')
        quit()
    
    if not os.path.isfile(args.TE_ref) :
        print('Error, cannot find TE reference file')
        quit()

    bam_name = args.bam
    bamfile = pysam.AlignmentFile(bam_name, "rb", reference_filename = args.ref)
    bam_header= bamfile.header
    
    bamfile.close( )

    try:
        sample_id=bam_header["RG"][0]["SM"]
    except:
        sample_id=bam_name.split("/")[-1].split(".")[0]

    samples=[sample_id]

    contigs = []
    contig_length = {}
    for contig in bam_header["SQ"]:
        contigs.append(contig["SN"])
        contig_length[ contig["SN"] ]=contig["LN"]

    vcf_header=tellr_vcf_header.main(bam_header, sample_id, version, contigs)
    TEs = tellr_call_vars.main(bamfile, samples, contigs) 

if __name__ == '__main__':
    main()
    