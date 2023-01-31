import sys
import pysam
import numpy as np

"""
fetch split reads
sort after split read size
get sequence
compare sequence to input file
"""

def get_region(bamfile, samples, chr, start, end, cand_dict, cand_id_dict, name_to_pos ):
    """
    extracts split reads from bam file to dictionaries
    """

    for read in bamfile.fetch(chr, start, end, until_eof=True):

        #print(chr, start, end)

        #only check reads that are mapped and with good quality
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        elif read.mapping_quality < 5:
            continue
        
        #get all split reads in region with support from more than X

        elif read.has_tag('SA') :
            #print(read)
            #print('qname', read.qname, 'pos', read.pos)

            cand_dict.append([chr, int(read.pos)])
            if chr not in cand_id_dict:
                cand_id_dict[chr] = {}
            cand_id_dict[chr][int(read.pos)] = read.qname
            name_to_pos[read.qname] = str(chr) + '-' + str(read.pos)

            #append chr and pos of supplementary read
            SR = read.get_tag('SA').rstrip(';').split(';')
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
    contigs = ['1', '4', '6']
    candidates = []
    candidates_toid = {}
    name_to_pos = {}
    for chr in contigs: 
        #if chr not in candidates:
            #candidates[chr] = []
        chr_len = contig_length[chr]
        for bin in range(1, chr_len-1000000, 1000000):
            get_region(bamfile, samples, chr, bin, bin+1000000, candidates, candidates_toid , name_to_pos) #returns candidates (dict with split read and its supplementaries) for clustering

    return candidates, candidates_toid, name_to_pos

