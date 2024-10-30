# sTELLeR

sTELLer; (s) Transposable ELements Long (e) Read, 
a tool for detection of non-reference transposable element using long-read data such as PacBio or ONT

![TELLR_flow](https://github.com/user-attachments/assets/ac9d5989-b56f-4356-8a91-7ab51ef6392d)

sTELLeR is a fast and CPU light tool for detection of transposable element insertion in long-read data. It is also possible to look for other insertions by giving the fasta sequence in --TE_fasta. 
sTELLeR gives output in VCF file, is haplotype-aware and can run on genome assemblies as well as on any species

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
require python3, pysam, samtools, minimap2 and miniasm

  python steller/steller.py [-h] [-R REF] -tf TE_FASTA -b BAM [-tr TE_REF] -s STYLE [-r SR] [-m MQ] [-k] [-o OUTPUT]

  arguments:

    -R REF, --ref     reference genome
    -tf, --TE_fasta   fasta file with sequence of elements to be detected
                          
    -b, --bam         bam file
                         
    -s , --style      ont or pb (clr) / hifi
                          
    -r , --sr         Minimum number of supporting split reads/insertions to call a variant (default 3)
    -m , --mq         Mapping quality (default 20)


  
  optional: 
  
    -k, --keep_intermediates Keep intermediate files
                        
    -o, --output      Output file name
                        

    -sr, --sr         Minimum number of supporting split reads to call a variant (default 3)

    -tr , --TE_ref    bed file with positions to avoid

NOTE: 
If using container, use python /sTELLeR/steller/steller.py 

Example to run on PB data: 

    python steller/steller.py --ref < genome ref_file > --TE_fasta < TE sequence fasta file > --bam < bamfile > --sr 4 --style pb -o < output prefix > -mr 80

Example to run on ONT data: 

    python steller/steller.py --ref < genome ref_file > --TE_fasta < TE sequence fasta file > --bam < bamfile > --sr 4 --style ont -o < output prefix > -mr 80

Example using test data located in: https://github.com/kristinebilgrav/sTELLeR_supplementary/
Download all files in repository and run: 
    python steller/steller.py --ref HG38_masked_noNN_chr22.fa --TE_fasta fasta/TEfastasequences.fa --bam testdata.bam --sr 4 --style pb -o testdata_res -mr 80
Results will be given in testdata_res_repeats.vcf and correct results can be verified by comparing to the provided file testbamTRUTH.txt - a total of 7 variants should be identified.  

# Output

The output is given in a VCF file (https://samtools.github.io/hts-specs/VCFv4.2.pdf), which is compatible to merge with VCF files from other variant callers. 
The variant is named INS:ME and in the info section the TE type along with start and end positions are given. 
In the format section (TE:TS:HT:GT) items such as TE type (as given in the fasta file), TS transposon start position, HT haplotag and GT genotype (as given in the bam file) are written out. 

If not specified otherwise, intermediate files containing analyzed reads (ID_chr__reads.txt) and candiate reads aligned to TEs (ID_chr_repeats.sam) will be removed. 

