import os
import gzip as gz
import pysam
import statistics

"""
maps sequence to TE and refines breakpoints
"""


def de_novo(chr, bam_name, repeat_fasta, sample, readfile, style):
    print('starting mapping')
    candidate_prefix = chr + '_' + sample + '_candidates'

    #extract reads of interest from bam to fasta
    bam_fasta = 'samtools view -h {} -N {} | samtools fasta -s stdout -n - > {}.fasta '.format(bam_name, str(readfile), candidate_prefix)
    print(bam_fasta)
    os.system(bam_fasta) 

    
    aligned_repeats = chr + '_' + sample + '_repeats.sam'

    #map fasta to TE fasta as short read 
    map_repeats = 'minimap2 -ax map-{} {}  {}.fasta >  {}'.format(style,repeat_fasta,candidate_prefix, aligned_repeats)
    print(map_repeats)
    os.system(map_repeats) 
    
    return aligned_repeats 



def check_cigar(flag, cigar, threshold):
    """
    reads cigar tag and counts bases to where
    the flag (type of var) is positioned
    """
    basesToAdd=[]

    num_flags = len([tup for tup in cigar if tup[0] == flag and tup[1] > threshold])

    if num_flags == 0:
        return False
    num_rounds = []

    add = 0
    if cigar[0][0] == flag:
        add += int(cigar[0][1])
    for c in cigar: 
        if c[0] != flag:
            add += c[1]

        else:   
            if c[1] > threshold:
                if c in num_rounds:
                    add += c[1]
                    continue
                num_rounds.append(c)
                basesToAdd.append(add)
                if len(num_rounds) == num_flags:
                    return basesToAdd
            else:
                continue

    return basesToAdd


def pos_toavoid(file): #change to pytabix
    """
    poisitions to avoid (known reference TEs) provided by user
    """

    PosToAvoid = {}
    for line in open(file):
        if line.startswith('#'):
            continue
        chr = line.split('\t')[0].strip('chr')
        if chr not in PosToAvoid:
            PosToAvoid[chr] = {}
        apos = line.split('\t')[1]
        bpos = line.split('\t')[2]

        PosToAvoid[chr][apos] = bpos
    return PosToAvoid


def extract_TEs(repeat_samfile):
    readToTEtype={} # Save readname and TE class
    readToTEPos={} # Save readname and position of identified TE

    avoid_flags = [2048, 2064]

    samfile = pysam.AlignmentFile(repeat_samfile, 'r')
    samfile_header = samfile.header
    TEs = {}
    for te in samfile_header["SQ"]:
        thete = te["SN"]
        length = te["LN"]
        TEs[thete] = length
    
    for line in samfile.fetch():

        # Only check approved reads
        flag = line.flag
        if flag in avoid_flags:
            continue

        read = line.qname # Readname
        r_start = line.pos # Get position of repeat in read - for crossreference with breakpoints
        repeat = line.reference_name # Repeat type
        if repeat not in TEs: # Confirm repeat in header
            continue

        # Save read to dictionaries
        if read not in readToTEPos:
            readToTEPos[read]={}
            readToTEtype[read]=[] # Can be multiple TEs in one read
        readToTEtype[read].append(repeat) # Add all types of TEs found in read

        if repeat not in readToTEPos[read]: # Save repeat start to repeat in read
            readToTEPos[read][repeat] = [] # Can be multiple TEs in one read
        readToTEPos[read][repeat].append(r_start)

     
    return readToTEtype, readToTEPos, TEs


def te_breakpoints(clusterToRead, readToTEtype, readToTEPos, clusterToPos , chr, TEs, haplodict, ReadStarts, ReadToVarPos):
    """
    Map confirmed TEs back to clusters
    Crosscheck that confirmed TE is same as clustered candidate
    //
    Go through each cluster 
    Check that reads in cluster are mapping to same TE
    Crosscheck that TE in read is same as clustered candidate
    """
    repeat_vars = []
    for c in clusterToRead:
        
        reads = clusterToRead[c]
        cluster_positions= list(clusterToPos[c]) # List of all positions in one cluster
        clusterconsensus = int(statistics.median(list(cluster_positions))) # Get middle position

        # Go through reads in cluster and check their TE match and haplotag
        thisclusterTEs={}
        for r in reads:
            if r not in readToTEtype:
                continue
           
            tes= readToTEtype[r]
            ht=haplodict[r]
            VarStart= ReadToVarPos[r]
            for v in VarStart:
                if v in cluster_positions:
                    for t in tes:
                        if t not in thisclusterTEs:
                            thisclusterTEs[str(t)]=[]
                        thisclusterTEs[str(t)].append(str(ht))

                    # Check that TE alignes with split reads and insertions
                    #pos_matches= [i for i in range(tepos-100, tepos+100) if i  in cluster_positions]
                    #if len(pos_matches) < 1:
                    #   continue
          
        for thete in thisclusterTEs:
            TEfrq=len(thisclusterTEs[thete])/len(reads)
            if TEfrq > 0.6:
                htags=list(set(thisclusterTEs[thete]))
                if len(htags) >1:
                    teht=','.join(htags)
                    tegt='1|1'
                else:
                    teht=htags[0]
                    if htags[0] == '1':
                        tegt='1|0'
                    elif htags[0] == '2':
                        tegt='0|1'
                    else:
                        teht='0'
                        tegt='0/0'

                lst = [thete, chr, str(clusterconsensus), str(clusterconsensus + int(TEs[thete])), teht, tegt ]
                if lst not in repeat_vars:
                    repeat_vars.append(lst)

            # Check that its not a reference TE
            # if chr in PosToAvoid:
                # start_positions = PosToAvoid[chr]
                # skip = [p for p in range(clusterconsensus-100, clusterconsensus+100) if p in start_positions ]#
                # if len(skip) > 0:
                    # continue
 
    return repeat_vars


def main(chr, bam_name, repeat_fasta, sample, readfile, clusterToPos, clusterToRead, refrepeat, style, haplotags, ReadStarts, ReadToVarPos):
    aligned = de_novo(chr, bam_name, repeat_fasta, sample, readfile, style)
    #    aligned = chr + '_' + sample + '_repeats.sam' 
    if refrepeat:
        avoid= pos_toavoid(refrepeat)

    tes=extract_TEs(aligned)
    readToTEtype =tes[0]
    readToTEPos = tes[1] 
    TEs = tes[2]
    variants = te_breakpoints(clusterToRead, readToTEtype, readToTEPos, clusterToPos, chr, TEs, haplotags, ReadStarts, ReadToVarPos)
    return variants

