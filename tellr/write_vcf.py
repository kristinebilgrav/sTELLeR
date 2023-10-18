import sys
import os

"""
write result to vcf
"""


def write_vcf(variants, filename, chr_lengths, sample):
    file =open(filename, 'w')
    vcfheader = ['##fileformat=VCFv4.1', '##source=TELLRv0',  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,REFME/NOVEL">' ]
    file.write('\n'.join(vcfheader) + '\n')

    for chr in chr_lengths:
        length = chr_lengths[chr]
        contigheader = '##contig=<ID={},length={}>'.format(chr,length)
        file.write(contigheader + '\n')
    
    formatheader =['##FORMAT=<ID=TS,Number=1,Type=String,Description="Transposon start">', '##FORMAT=<ID=TE,Number=1,Type=String,Description="Transposable element">', 
                   '##FORMAT=<ID=HT,Number=1,Type=String,Description="Haplotag">', '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    file.write('\n'.join(formatheader) + '\n')

    vcfheader2 = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample]
    file.write('\t'.join(vcfheader2) + '\n')

    for entry in variants:
        for subentry in entry:
            te = subentry[0]
            chr = subentry[1]
            start = subentry[2]
            end = subentry[3]
            ht = subentry[4]
            gt=subentry[5]
            meinfo = ','.join(['MEINFO=' + te, str(start), str(end)])
            svtype = 'SVTPE=INS' 
            infoline = ';'.join([meinfo, svtype])
            sampleformat = ':'.join([te,str(start),ht,gt ])
            outline = [chr, str(start), '.', '.', '<INS:ME>', '.', '.', infoline, 'TE:TS:HT:GT', sampleformat ]
            file.write('\t'.join(outline) + '\n') 
            
    return file

def del_files(sample,chrs ):
    os.system('rm {}_reads.txt'.format(sample))
    for chr in chrs: 
        candidate_prefix = chr + '_' + sample + '_candidates'
        os.system('rm {}.fasta'.format(candidate_prefix))
        aligned_repeats = chr + '_' + sample + '_repeats.sam'
        os.system('rm {}'.format(aligned_repeats))


def main(variants, chr_lengths, sample, chrs):
    var_out = sample + '_repeats.vcf'
    write_vcf(variants, var_out,chr_lengths ,sample)
    del_files(sample, chrs)