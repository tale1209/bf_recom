#!/bin/sh

wd="path/to/wd" && cd $wd
bubble_matrix="path/to/bubble_matrix"     #bubble_matrix file, single chr     #delimiter:'\t'
output_file="path/to/output"

missing=0.2  
binomial_p=0.05

min_sample_number=$(awk 'NR!=1{print NF-1}' $bubble_matrix | sort | uniq | awk '{print $1*(1-'$missing')}')

if [ -e $output_file ];then rm $output_file;fi && awk 'NR!=1' $bubble_matrix | while read line;do
  ps_number=$(echo $line | awk '{c=0;for(i=2;i<=NF;i++)if($i!="NA")c++;print c}')
  if echo "$ps_number >= $min_sample_number" | bc -l | grep -q 1;then      #filter1:missing rate
    p_number=$(echo $line | awk '{c=0;for(i=2;i<=NF;i++)if($i=="0")c++;print c}')        # 0-primary 1-secondary
    binomial=$(python ${wd}/binomial.py -k $p_number -n $ps_number | head -n 4 | tail -n 1 | sed 's/P-value://g')
    if echo "$binomial >= $binomial_p" | bc -l | grep -q 1;then       #filter2:binomial test
      echo $line | sed 's/ /\t/g' >> $output_file
    fi
  fi
done

cat <(head -n 1 $bubble_matrix) $output_file > tmp && mv -f tmp $output_file
