#!/bin/sh

wd="path/to/wd" && cd $wd
ref_chr="path/to/ref_fna"
hap_chr="path/to/hap_fna"     #according to hapi phasing

minimap2 -a -x asm20 --cs -r2k -t 5 $ref_chr $hap_chr > ref_vs_hap.sam
samtools sort -m4G -@4 -o ref_vs_hap.sorted.bam ref_vs_hap.sam && samtools index ref_vs_hap.sorted.bam
svim-asm haploid ./ ref_vs_hap.sorted.bam $ref_chr --query_names
