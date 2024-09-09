from multiprocessing import Pool
import argparse
import pathlib
import os
import pysam

import variants as steller_call_vars
import cluster as steller_cluster_vars
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
        prog = 'sTELLeR', 
        description='calls non-reference transposable elements given in a long-read pacbio or ont bam file'
    )

    parser.add_argument('-R', '--ref', help='reference genome')
    parser.add_argument('-tf', '--TE_fasta', help='fasta file with elements to be detected', required=True, type=pathlib.Path)
    parser.add_argument('-b', '--bam', help='bam file', required= True, type=pathlib.Path)
    parser.add_argument('-ta', '--TE_avoid', help='bed file with positions to avoid', required= False, type=pathlib.Path)
    parser.add_argument('-s', '--style', help='ont or pb', required= True, type=pathlib.Path)
    parser.add_argument('-r', '--sr', help='Minimum number of supporting split reads/insertions to call a variant (default 3)', required= False, default = 3)
    parser.add_argument('-mr', '--maxreads', help='Maximum number of supporting split reads/insertions to call a variant (default 100)', required= False, default = 100)
    parser.add_argument('-m', '--mq', help='Mapping quality (default 20)', required= False, default = 20)
    parser.add_argument('-k', '--keep_intermediates', help='Keep intermediate files', required= False, action="store_false")
    parser.add_argument('-o', '--output', help='Output file name', required= False)

    args = parser.parse_args()

    """
    check files
    """
    if not os.path.isfile(args.bam) :
        print('Error, cannot find bam file')
        quit()
    
    elif not os.path.isfile(args.TE_fasta) :
        print('Error, cannot find TE fasta file')
        quit()

    else:
        parser.print_help()

    """
    define parameters needed in modules
    """
    repeat_fasta = args.TE_fasta
    sr = int(args.sr)
    mr = int(args.maxreads)
    mapping_quality = int(args.mq)
    style = args.style
    bam_name = str(args.bam)
    bamfile = pysam.AlignmentFile(bam_name, "rb", reference_filename = args.ref)
    bam_header= bamfile.header

    if args.output:
        sample_id=args.output
    else:
        sample_id=bam_name.split("/")[-1].split(".bam")[0]
    sample=sample_id

    if args.TE_avoid is None :
        repeatsToAvoid = False
    else:
        repeatsToAvoid = args.TE_avoid
    

    delete= args.keep_intermediates

    
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
        candidates = steller_call_vars.main(chr, bamfile, chr_length, mapping_quality) 
        chr_candidates = candidates[0]
        chr_candidatereads=candidates[1]
        readinfo= candidates[2]


        print('clustering')
        clustered = steller_cluster_vars.main(chr, chr_candidates, chr_candidatereads, readinfo, sample, sr, mr )
        if clustered == False: 
            continue

        readsfile=clustered[1]
        clusterToRead=clustered[2]

        print('calling')
        calls = assembleandcall.main(chr, bam_name, repeat_fasta, sample, readsfile, clusterToRead, repeatsToAvoid, style, readinfo, sr) 
        repeatvariants.append(list(calls))

    print('writing to file')
    write_calls.main(repeatvariants, chr_length, sample, chrs, delete)



if __name__ == '__main__':
    main()
    