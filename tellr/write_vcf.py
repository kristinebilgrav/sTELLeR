import sys


def write_vcf(variants, filename):
    file =open(filename, 'w')
    for entry in variants:
        for subentry in entry:
            file.write('\t'.join(subentry) + '\n')
    return file


def main(variants, chr, sample):
    var_out = sample + '_repeats.vcf'
    write_vcf(variants, var_out)