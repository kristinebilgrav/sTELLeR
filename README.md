# sTELLeR

**sTELLeR; (s) Transposable ELements Long (e) Read**

a tool for detection of non-reference transposable element using long-read data such as PacBio or ONT

![TELLR_flow](https://github.com/user-attachments/assets/ac9d5989-b56f-4356-8a91-7ab51ef6392d)

**sTELLeR** is a fast and CPU light tool designed and tested for detection of transposable element insertion in long-read data, however it supports identification of any type of insertions. 
sTELLeR gives output in VCF, is haplotype-aware and can run on genome assemblies as well as on any species

Example scripts, bam test file and TE fasta sequences are available at https://github.com/kristinebilgrav/sTELLeR_supplementary/

# Install

*Option 1:*

git clone https://github.com/kristinebilgrav/sTELLeR.git

Install dependencies: 
- minimap2 
- pysam
- samtools 

All available through bioconda. 

*Option 2:* 

Dowload docker container containing sTELLeR available at:
https://hub.docker.com/r/kristinebilgrav/steller

run using docker or singularity with:

    python /sTELLeR/steller/steller.py 

# Run

sTELLeR needs a bam or cram file, the reference used to align it, and a fasta file with insertions sequences to look for (--TE_fasta). This allows for detection of any type of insertions, not only TE insertions.
Example TE sequences are given in the test dataset [HERE](https://github.com/kristinebilgrav/sTELLeR_supplementary/) \
If the input is phased, the output VCF will contain information regarding which haplotype(s) the TE is on and its genotype. If not phased the haplotype will be 0 and the genotype will indicate if TE is heterozygous (1/0) or homozygous (1/1).

require python3, pysam, samtools, minimap2 

    python steller/steller.py [-h] [-R REF] -tf TE_FASTA -b BAM [-tr TE_REF] -s STYLE [-r SR] [-m MQ] [-k] [-o OUTPUT]

  arguments:

    -R REF, --ref     reference genome
    -tf, --TE_fasta   fasta file with sequence of elements to be detected
                          
    -b, --bam         bam file
                         
    -s , --style      ont or pb (clr) / hifi
                          
  
  optional: 
  
    -k, --keep_intermediates Keep intermediate files
                        
    -o, --output      Output file name
    
    -r , --sr         Minimum number of supporting split reads/insertions to call a variant (default 3)
    
    -m , --mq         Mapping quality (default 20)
                        
    -mr, --maxreads   Maximum number of supporting split reads/insertions to call a variant (default 100)


The intermediate files (--keep_intermediates) include per-chromosome fasta files with candidate insertions as well as SAM files with aligned insertions. \
Supporting reads (--sr) of a insertion depends on the coverage of the sample and the desired sensitivity.\
Mapping quality (--mq) is a threshold of accepted TE calls to the specific TE sequence. \
Maximum number of reads (--maxreads) can be adjusted to avoid calls in outlying regions. 



**NOTE:**
If using container, the executable is:

    python /sTELLeR/steller/steller.py 


### Examples 

Example to run on PB data: 

    python steller/steller.py --ref < genome ref_file > --TE_fasta < TE sequence fasta file > --bam < bamfile > --sr 4 --style pb -o < output prefix > -mr 80

Example to run on ONT data: 

    python steller/steller.py --ref < genome ref_file > --TE_fasta < TE sequence fasta file > --bam < bamfile > --sr 4 --style ont -o < output prefix > -mr 80


Example using test data located in: https://github.com/kristinebilgrav/sTELLeR_supplementary/

Download all files in repository and run: 

    python steller/steller.py --ref HG38_chr22.fa --TE_fasta fasta/TEfastasequences.fa --bam testdata.bam --sr 4 --style pb -o testdata_res -mr 80

Results will be given in testdata_res_repeats.vcf and correct results can be verified by comparing to the provided file testbamTRUTH.txt - a total of 7 variants should be identified.  

# Output

The output is given in a [VCF file](https://samtools.github.io/hts-specs/VCFv4.2.pdf), which is compatible to merge with VCF files from other variant callers. \
The variant is named INS:ME and in the info section the TE type along with start and end positions are given. \
In the format section (TE:TS:HT:GT) items such as TE type (as given in the fasta file), TS transposon start position, HT haplotag (0 if not phased) and GT genotype (1/0 or 1/1 for unphased and 1|0, 0|1 or 1|1 for phased) are written out. 

If not specified otherwise, intermediate files containing fasta sequences of candidate insertions (ID_chr_reads.fasta) and candiate reads aligned to TEs (ID_chr_repeats.sam) will be removed. 

