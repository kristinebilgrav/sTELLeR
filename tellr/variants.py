import sys
import pysam
import numpy as np

"""
fetch split reads
sort after split read size
get sequence
compare sequence to input file
"""

def check_cigar(flag, cigar, read, threshold):
    basesToAdd=0

    #num_flags = len([tup for tup in cigar if tup[0] == flag and tup[1] > threshold])
    #print(num_flags)


    


    for c in cigar: 
        if c[0] != flag:
            basesToAdd += c[1]

            continue
        else:   
            if c[1] > threshold:
                break
            else:
                continue

    return basesToAdd

def get_region(bamfile, samples, chr, start, end, cand_dict, cand_id_dict, name_to_pos ):
    """
    extracts split reads from bam file to dictionaries
    """

    for read in bamfile.fetch(chr, start, end, until_eof=True):
        chr = chr.strip('chr')
        if chr not in cand_id_dict:
            cand_id_dict[chr] = {}

        #print(chr, start, end)

        #only check reads that are mapped and with good quality
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        elif read.mapping_quality < 5:
            continue
        
        #get all split reads in region with support from more than X

        elif read.has_tag('SA') :
            #explain_cigar = 0: matching, 1:insertion, 2:deletion, 3:ref_skip, 4:soft_clipped, 5:hard_clipped, 6:...
            

            #find split read position
            #reverse complement - oposite of others
            cigg = read.cigar
            #print(cigg)
            read_start_pos = read.pos
            soft_clipping = check_cigar(4, cigg, read, 100)
            #print('soft')


            insertions = check_cigar(1, cigg, read, 100)
            #print('insertions')

            variants = [soft_clipping, insertions]
            for v in variants:

                splitpos = read_start_pos + v

                if read.is_reverse:
                    splitpos = read_start_pos - v

                #print(read_start_pos, splitpos)
                cand_dict.append([chr, splitpos])
                if splitpos not in cand_id_dict[chr]:
                    cand_id_dict[chr][splitpos] = []

                cand_id_dict[chr][splitpos].append(read.qname)
                name_to_pos[read.qname] = str(chr) + '-' + str(splitpos)

            #append chr and pos of supplementary read
            #SR = read.get_tag('SA').rstrip(';').split(';')
            #for sa in SR:

                #if sa_chr not in cand_id_dict :
                #    cand_id_dict[sa_chr] = {}
                #sa_pos = int(sa.split(',')[1])
                #cand_dict.append([sa_chr, sa_pos])
                #cand_id_dict[sa_chr][sa_pos] = 

        else:
            continue

    return cand_dict, cand_id_dict, name_to_pos
 



def main(bamfile, samples, contigs, contig_length, sr):
    #TEs = []
    contigs = ['chr6']
    candidates = [] #[chr, pos]
    candidates_toid = {} #chr:{pos:[readname(s)]}
    name_to_pos = {} #readname:chr-pos
    for chr in contigs: 
        #if chr not in candidates:
            #candidates[chr] = []
        chr_len = contig_length[chr.strip('chr')]
        for bin in range(1, chr_len-1000000, 1000000):
            get_region(bamfile, samples, chr, bin, bin+1000000, candidates, candidates_toid , name_to_pos) #returns candidates (dict with split read and its supplementaries) for clustering

    return candidates, candidates_toid, name_to_pos

