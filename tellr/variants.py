import sys
import pysam
import numpy as np


"""
fetch split reads
sort after split read size
get sequence
compare sequence to input file
"""

def check_cigar(flag, cigar, threshold):
    """
    Reads cigar tag and counts bases to where
    the flag (type of var) is positioned
    threshold: needs to be at least threshold amount of soft clipped / inserted bases
    """
    RefbasesToAdd={}
    Seqbases={}

    num_flags = len([tup for tup in cigar if tup[0] == flag and tup[1] > threshold]) # For each pair (flag,bases) in cigar checks 0: flag and 1:bases to see if desired flag and >threshold

    if num_flags == 0:
        return False
    
    num_rounds =[] # Length of this to be same as length of num_flags
    avoid = [1,5]
    seqavoid=[2, 3, 5]
    add = 0     
    seqadd=0
    for c in cigar:
        if c[0] != flag:
            if c[0] not in avoid: # Dont add bases for insertions and hard clipped to find RefPos
                add += c[1]
            if c[0] not in seqavoid: # Dont add bases for hard clipped and deletions, do for insertions
                seqadd += c[1] 

        else: # Matches flag
            if int(c[1]) > threshold: # Threshold                 
                num_rounds.append(c)
            
                if add not in RefbasesToAdd:

                    RefbasesToAdd[add]=c[1]
                    Seqbases[add]=seqadd


                if len(num_rounds) == num_flags:
                    return RefbasesToAdd, Seqbases
            else:
                if c[0] not in avoid: # Dont want to add bases for insertions or hard clipped reads as output is reference base
                    add += c[1]
                if c[0] not in seqavoid:
                    seqadd += c[1]
                
                continue

    return RefbasesToAdd, Seqbases

def get_region(bamfile, chr, start, end, cand_pos, cand_reads, readinfo, mapping_quality):
    """
    Runs through bam file and extracts positions of split reads and insertions in primary reads with mapping quality > user implied quality
    """
    for read in bamfile.fetch(chr, start, end, until_eof=True):

        # Only check reads that are mapped and with good quality
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        elif read.mapping_quality < mapping_quality: # Mapping quality according to user set quality
            continue

        # Retrieve needed variables
        cigg = read.cigartuples # Cigar: [ (1,24),(7,133) ]
        read_start_pos = read.reference_start # Read Start
        readname=read.qname

        # Explain_cigar = 0: matching, 1:insertion, 2:deletion, 3:ref_skip, 4:soft_clipped, 5:hard_clipped, 6:padding, 7:sequence match, 8:sequence mismatch
        # Check for soft clipped 
        soft_clipping = check_cigar(4, cigg, 100) #list w bases to add for each soft clip OR false
        # Check for insertions
        insertions = check_cigar(1, cigg, 10) # list w bases to add for each ins OR false
 
        # Either could be a TE 
        if soft_clipping or insertions :
      
            # Find split read position
            # Reverse complement - opposite of others?
            
            variants = [soft_clipping, insertions]
            for v in variants:
                if v == False:
                    continue
                for add in v[0]:


                    varpos = read_start_pos + add
                    varpos_len= v[0][add]
                    seqs=v[1][add]
                    if seqs > 100:
                        seqs = seqs-100
                    seqe= v[1][add] + varpos_len
                    if seqe+100 < read.query_length:
                        seqe= seqe+100
                    sequence= read.query_sequence[seqs:seqe]
                    orientation = '+'

                    #if len(sequence)<1:
                    #    print(read)
                    #   print(seqs, seqe, read.query_length)
                        #print(sequence)
                        

                    # Add position to list for clustering
                    cand_pos.append(varpos)

                    # Add readname to list/dict for later retrieval
                    cand_reads.append(readname)

                    # Retrive HP tag for phasing/genotyping
                    try:
                        haplotype = read.get_tag('HP') 
                    except:
                        haplotype = 0
        
                    # Add read info to dictionary
                    readdictname= readname+'_'+str(varpos)
                    # readstart, variantpos, length, orientation, haplotag, sequence
                    infotup= (read_start_pos,varpos, varpos_len, orientation,haplotype, sequence)
                    if readdictname not in readinfo:
                        readinfo[readdictname]=infotup

            
        else:
            continue

    return cand_pos, cand_reads, readinfo
 


def main(chr, bamfile, contig_length, mapping_quality):
    candidates = [] # [pos] return for each chromosome
    candidate_readIDs=[]
    readinfo = {} # Dict with key readname_varpos:(readstart, variantpos, length, orientation, haplotag, sequence)

    chr_len = contig_length[chr]
    for bin in range(1, chr_len-1000000, 1000000):
        get_region(bamfile, chr, bin, bin+1000000, candidates,candidate_readIDs, readinfo, mapping_quality) # Returns candidates (dict with split read and its supplementaries) for clustering
    

    return candidates, candidate_readIDs, readinfo
