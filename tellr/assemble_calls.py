import os
import gzip as gz

def de_novo(bam_name, repeat_fasta, sample, readfile):
    candidate_prefix = sample + '_candidates'
    candidate_bam = 'samtools view -h {} -N {} -b -o {}.bam'.format(bam_name, str(readfile), candidate_prefix)
    print(candidate_bam)
    os.system(candidate_bam) #save all candidate reads to bam

    bam_fasta = 'samtools fasta  {}.bam > {}.fasta '.format(candidate_prefix, candidate_prefix)
    print(bam_fasta)
    os.system(bam_fasta) #bam to fasta

    paf_file = sample + '_candidate_overlap.paf.gz'
    map_self = 'minimap2 -x ava-pb {}.fasta {}.fasta| gzip -1 > {}'.format(candidate_prefix, candidate_prefix, paf_file)
    print(map_self)
    os.system(map_self) #map reads against themself, creating read overlap 

    gfa_file = sample + '_candidate_overlap.gfa'
    denovo_assembly = 'miniasm -m 300 -f {}.fasta {} > {}'.format(candidate_prefix, paf_file, gfa_file)
    print(denovo_assembly)
    os.system(denovo_assembly) #de novo assembly

    denovo_fasta = sample + '_denovo.fasta'
    gfa_fasta= """ awk '/^S/{{print ">"$2"\\n"$3}}' {0} > {1} """.format(gfa_file, denovo_fasta)
    print(gfa_fasta)
    os.system(gfa_fasta) #gfa to fasta

    repeat_paffile = sample + '_candidate_repeat_overlap.paf.gz'
    #map_repeats = 'minimap2 -x ava-pb -c --cs {}  {} | gzip -1 >  {}'.format(denovo_fasta , repeat_fasta, repeat_paffile)
    map_repeats = 'minimap2 -ax sr {}  {} >  {}'.format(denovo_fasta , repeat_fasta, sample + '_candidates.bam')
    print(map_repeats)
    os.system(map_repeats) #map de novo fasta to TE fasta
    

    return gfa_file, repeat_paffile


def breakpoints(gfa_file, repeat_paffile, readNamedict, sample, clusterPosToName, PosToAvoid):
    
    #extracting de novo contig names and original read names
    denovoToReadName = {}
    for overlap in open(gfa_file):
        if not overlap.startswith('a') :
            continue

        denovo= overlap.split('\t')[1]
        reads = overlap.split('\t')[3].split(':')[0]
        if denovo not in denovoToReadName:
            denovoToReadName[denovo] = []
        denovoToReadName[denovo].append(reads)

        # paftools.js

    #repeats and contig matching more than 90% of repeat length
    #saved to dict denovo:repeat
    if 'gz' in repeat_paffile:
        opener = gz.open
    else:
        opener =open

    denovoToRepeat = {}
    var_out = open(sample + '_repeats.vcf', 'w')
    for line in opener(repeat_paffile, 'rt'):

        r_len = int(line.split('\t')[1])
        match = int(line.split('\t')[9])
        perMatch = match/r_len
        if perMatch < 0.7:
            continue

        repeat = line.split('\t')[0]
        #print(repeat)
        denovo_contig = line.split('\t')[5]
        if denovo_contig not in denovoToRepeat:
            denovoToRepeat[denovo_contig] = []
        denovoToRepeat[denovo_contig].append(repeat)

        reads = denovoToReadName[denovo_contig]
        #print(reads)

        chr = ''
        break_pos = 0
        for r in reads:
            chr = clusterPosToName[r].split('-')[0]
            pos = int(clusterPosToName[r].split('-')[1])
            break_pos = pos
            #print(pos)

        #check that its not a reference TE
        if chr in PosToAvoid:
            for spos in PosToAvoid[chr]:
                if break_pos in range(spos-100, PosToAvoid[chr][pos]+100) : #make into TE specific?
                    continue

        
        var_out.write('\t'.join([repeat, chr, str(break_pos), str(perMatch)]) + '\n') #add to set first and write second -should not be more elements than clusters. 
  

    return denovoToRepeat, denovoToReadName

def pos_toavoid(file):

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

def main(bam_name, repeat_fasta, sample, readfile, readName_dict, clusterPosToName, refrepeat):
    #aligned = de_novo(bam_name, repeat_fasta, sample, readfile)
    aligned = [sample + '_candidate_overlap.gfa',  sample + '_candidate_repeat_overlap.paf.gz'] 
    avoid= pos_toavoid(refrepeat)
    breakpoints(aligned[0], aligned[1], readName_dict, sample, clusterPosToName, avoid)
    #write vcf?

