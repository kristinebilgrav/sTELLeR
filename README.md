# TELLR

non-reference transposable element calling on long-read data 

# Run
require python3, pysam, samtools, minimap2 and miniasm

  python tellr/tellr.py [-h] [-R REF] -tf TE_FASTA -b BAM [-tr TE_REF] -s STYLE [-r SR] [-m MQ] [-k] [-o OUTPUT]

  arguments:

    -R REF, --ref     reference genome
    -tf, --TE_fasta   fasta file with sequence of elements to be detected
                          
    -b, --bam         bam file
    -tr, --TE_ref     bed file with positions to avoid
                         
    -s , --style      ont or pb 
                          
    -r , --sr         Minimum number of supporting split reads/insertions to call a variant (default 3)
    -m , --mq         Mapping quality (default 20)


  
  optional: 
    -k, --keep_intermediates Keep intermediate files
                        
    -o, --output      Output file name
                        

    -sr, --sr         Minimum number of supporting split reads to call a variant (default 3)

    -tr , --TE_ref    bed file with positions to avoid
                          

