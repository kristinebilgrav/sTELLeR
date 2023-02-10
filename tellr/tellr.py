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


import vcf_header as tellr_vcf_header 
import variants as tellr_call_vars
import cluster as tellr_cluster_vars
import assemble_calls as assembleandcall

def main():
    version = "0.0.0"
    parser= argparse.ArgumentParser(
        prog = 'tellr', 
        description='calls non-reference transposable elements given in xx file in long-read pacbio or ont bam files'
    )

    parser.add_argument('-ref', '--ref', help='reference genome')
    parser.add_argument('-TE_fasta', '--TE_fasta', help='fasta file with elements to be detected', required=True, type=pathlib.Path)
    parser.add_argument('-bam', '--bam', help='bam file', required= True, type=pathlib.Path)
    parser.add_argument('-sr', '--sr', help='Minimum number of supporting split reads to call a variant (default 3)', required= False, default = 3)
    parser.add_argument('-TE_ref', '--TE_ref', help='bed file with positions to avoid', required= True, type=pathlib.Path)

    args = parser.parse_args()

    if not os.path.isfile(args.bam) :
        print('Error, cannot find bam file')
        quit()
    
    elif not os.path.isfile(args.TE_ref) :
        print('Error, cannot find TE reference file')
        quit()

    else:
        parser.print_help()

    repeat_fasta = args.TE_fasta
    repeatsToAvoid = args.TE_ref
    bam_name = str(args.bam)
    bamfile = pysam.AlignmentFile(bam_name, "rb", reference_filename = args.ref)
    bam_header= bamfile.header
    print(bam_header)
    
    #bamfile.close()


    sample_id=bam_name.split("/")[-1].split(".")[0]

    sample=sample_id

    chrs = []
    chr_length = {}
    for contig in bam_header["SQ"]:
        thischr = contig["SN"].strip('chr')
        chrs.append(thischr)
        chr_length[thischr]=contig["LN"]


    #vcf_header=tellr_vcf_header.main(bam_header, sample_id, version, contigs)
    candidates = tellr_call_vars.main(bamfile, sample, chrs, chr_length, args.sr) 
    clustered = tellr_cluster_vars.main(candidates[0], candidates[1], bamfile, sample, bam_name, repeat_fasta )
    calls = assembleandcall.main(bam_name, repeat_fasta, sample, clustered[1], candidates[2], clustered[2], repeatsToAvoid) 

if __name__ == '__main__':
    main()
    