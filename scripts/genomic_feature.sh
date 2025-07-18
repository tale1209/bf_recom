#!/bin/sh

wd="path/to/wd" && cd $wd
ref="path/to/chromosome-level-genome"              #reference
crossover="path/to/co_result"        #merged crossover file generated by hapi, replace the 1st column with the corresponding chromosome
coldspot="path/to/coldspot"         #coldspot file
anno="path/to/genomic_feature_files"             #annotation files(Gene/TE...)

bin_size=5000000
step=`expr $bin_size / 5`

python ${wd}/correlation.py pearson \
  <(bedtools coverage \
    -a <(awk -v OFS='\t' 'NR!=1{print $1,$5,$5}' $crossover | bedtools coverage \
          -a <(bedtools makewindows -g "$ref".fai -w $bin_size -s $step | sort -k1,1V -k2,2n -k3,3n | bedtools subtract -a - -b $coldspot) \
          -b - -counts) \
    -b <(awk -v OFS='\t' '{print $1,int(($3+$2)/2),int(($3+$2)/2)}' $anno) -counts | sed '1ichr\tstart\tend\tvalue1\tvalue2') \
value2 value1


##GC content
python ${wd}/correlation.py pearson \
  <(bedtools coverage \
    -a <(bedtools makewindows -g "$ref".fai -w $bin_size -s $step | sort -k1,1V -k2,2n -k3,3n | bedtools subtract -a - -b $coldspot \
         | bedtools getfasta -fi $ref -bed - | seqkit fx2tab -l -g -n -i -H | grep -v "#" \
         | awk -v OFS='\t' -F '[\t: -]' '{print $1,$2,$3,$5}') \
    -b <(awk -v OFS='\t' 'NR!=1{print $1,$5,$5}' $crossover) -counts \
    | awk -v OFS='\t' '{print $1,$2,$3,$5,$4}' | sed '1ichr\tstart\tend\tvalue1\tvalue2') \
value2 value1
