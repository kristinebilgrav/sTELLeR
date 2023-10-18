import sys
import argparse
import pathlib
import os
import pysam
import variants as tellr_call_vars
import cluster as tellr_cluster_vars
import assemble_calls as assembleandcall
import write_vcf as write_calls


"""
import py script modules
define main
parsers
setups/wrapper
"""

def main():
    version = "0.0.0"
    parser= argparse.ArgumentParser(
        prog = 'tellr', 
        description='calls non-reference transposable elements given in a long-read pacbio or ont bam file'
    )

    parser.add_argument('-R', '--ref', help='reference genome')
    parser.add_argument('-tf', '--TE_fasta', help='fasta file with elements to be detected', required=True, type=pathlib.Path)
    parser.add_argument('-b', '--bam', help='bam file', required= True, type=pathlib.Path)
    parser.add_argument('-tr', '--TE_ref', help='bed file with positions to avoid', required= True, type=pathlib.Path)
    parser.add_argument('-s', '--style', help='ont or pb', required= True, type=pathlib.Path)
    parser.add_argument('-r', '--sr', help='Minimum number of supporting split reads/insertions to call a variant (default 3)', required= False, default = 3)
    parser.add_argument('-m', '--mq', help='Mapping quality (default 20)', required= False, default = 20)

    args = parser.parse_args()

    """
    check files
    """
    if not os.path.isfile(args.bam) :
        print('Error, cannot find bam file')
        quit()
    
    elif not os.path.isfile(args.TE_ref) :
        print('Error, cannot find TE reference file')
        quit()

    else:
        parser.print_help()

    """
    define parameters needed in modules
    """
    repeat_fasta = args.TE_fasta
    sr = int(args.sr)
    mapping_quality = int(args.mq)
    style = args.style
    repeatsToAvoid = args.TE_ref
    bam_name = str(args.bam)
    bamfile = pysam.AlignmentFile(bam_name, "rb", reference_filename = args.ref)
    bam_header= bamfile.header
    sample_id=bam_name.split("/")[-1].split(".")[0]
    sample=sample_id
    
    #bamfile.close()

    # Note down chromosomes and their lenght from the bam file
    chrs = []
    chr_length = {}
    for contig in bam_header["SQ"]:
        thischr = contig["SN"]
        chrs.append(thischr)
        chr_length[thischr]=contig["LN"]
    

    """
    for each chromosome; find candidates, cluster, align to TE and collect positions
    write to vcf when performed on allp
    """
    repeatvariants =[]
    for chr in chrs:
        print('finding candidates on', chr)        
        candidates = tellr_call_vars.main(chr, bamfile, chr_length, mapping_quality) 
        chr_candidates = candidates[0]
        candidates_toid = candidates[1]
        ReadStarts=candidates[2]
        HaploTags = candidates[3]
        ReadToVarPos = candidates[4]

        print('clustering')
        clustered = tellr_cluster_vars.main(chr, chr_candidates, candidates_toid, bamfile, sample, bam_name, sr, repeat_fasta )
        if clustered == False: 
            continue

        readtxtfile=clustered[1]
        clusterToPos=clustered[3]
        clusterToRead=clustered[2]

        print('calling')
        calls = assembleandcall.main(chr, bam_name, repeat_fasta, sample, readtxtfile, clusterToPos, clusterToRead, repeatsToAvoid, style, HaploTags, ReadStarts, ReadToVarPos) 
        repeatvariants.append(list(calls))

    print('writing to file')
    write_calls.main(repeatvariants, chr_length, sample, chrs)


if __name__ == '__main__':
    main()
    