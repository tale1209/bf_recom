#!/bin/sh

wd="path/to/wd" && cd $wd
ref_pm="path/to/parental_fna"            # parental genome
ref_p="path/to/paternal_fna"             # paternal genome
ref_m="path/to/maternal_fna"             # maternal genome

# identify difference in echo bubble pair
parents="P"           # P or M
bubble="1"            # all bubbles needed
awk '$1~"primary" && $1~"bubble'$bubble'_" && $1~"Bf_'$parents'"' $ref_pm | sed 's/>//g' | seqkit grep -f - $ref_pm > "$parents"_bubble"$bubble"_pri.fa
awk '$1~"secondary" && $1~"bubble'$bubble'_" && $1~"Bf_'$parents'"' $ref_pm | sed 's/>//g' | seqkit grep -f - $ref_pm > "$parents"_bubble"$bubble"_sec.fa
nucmer --mum -t 3 --prefix "$parents"_bubble"$bubble" "$parents"_bubble"$bubble"_pri.fa "$parents"_bubble"$bubble"_sec.fa
dnadiff -p "$parents"_bubble"$bubble" -d "$parents"_bubble"$bubble".delta       # primary VS secondary 

bin_size=100
step=`expr $bin_size / 5`
cat "$parents"_bubble*.snps | awk '{print $9,$11,$1,$2,$3,$4,$12,$10}' > pvss_diff_"$parents"
# 1-5 SNPs per 100bp
awk 'NR==FNR{a[$2]=1;next}$1 in a{print $0}' pvss_diff_"$parents" <(bedtools makewindows -g "$ref_pm".fai -w $bin_size -s $step) \
| bedtools coverage -a - -b <(awk -v OFS='\t' '{print $2,$3,$3}' pvss_diff_"$parents") -counts \
| awk -v OFS='\t' '$4>=1 && $4<=5{print $1,$2,$3}' | bedtools merge \
| bedtools intersect -a <(awk -v OFS='\t' '{print $2,$3,$3}' pvss_diff_"$parents") -b - > "$parents"_marker
# remove indels
awk '$4!="." && $4!="N" && $5!="." && $5!="N"{print $0}' pvss_diff_"$parents" \
| awk 'NR==FNR{a[$2"_"$3]=1;next}$1"_"$2 in a{print $0}' - "$parents"_marker > tmp && mv -f tmp "$parents"_marker
# merge primary and secondary
awk -v OFS='\t' 'NR==FNR{a[$1"_"$2]=1;next}$2"_"$3 in a{print $7,$6,$6}' "$parents"_marker pvss_diff_"$parents" \
| cat "$parents"_marker - > tmp && mv -f tmp "$parents"_marker
#excluded parental homologous markers
nucmer nucmer --mum -t 3 --prefix PvsM $ref_p $ref_m && dnadiff -p PvsM -d PvsM.delta  
if [ $parents == "P" ];then
awk '{print $11"_"$1"_"$2"_"$3}' PvsM.snps \
| awk -v OFS='\t' 'NR==FNR{a[$1]=1;next}!($1"_"$2"_"$3"_"$4 in a){print $1,$2,$2}' - \
<(awk 'NR==FNR{a[$2"_"$3]=$4 FS $5;b[$7"_"$6]=$5 FS $4;next}$1"_"$2 in a{print $1,$2,a[$1"_"$2]}$1"_"$2 in b{print $1,$2,b[$1"_"$2]}' pvss_diff_"$parents" "$parents"_marker) \
> tmp && mv -f tmp "$parents"_marker
elif [ $parents == "M" ];then
awk '{print $12"_"$4"_"$3"_"$2}' PvsM.snps \
| awk -v OFS='\t' 'NR==FNR{a[$1]=1;next}!($1"_"$2"_"$3"_"$4 in a){print $1,$2,$2}' - \
<(awk 'NR==FNR{a[$2"_"$3]=$4 FS $5;b[$7"_"$6]=$5 FS $4;next}$1"_"$2 in a{print $1,$2,a[$1"_"$2]}$1"_"$2 in b{print $1,$2,b[$1"_"$2]}' pvss_diff_"$parents" "$parents"_marker) \
> tmp && mv -f tmp "$parents"_marker
else
echo "wrong mode"
fi
#rm contigs with <3 markers
awk '{print $1}' "$parents"_marker | sort | uniq -c | awk '$1<3{print $2}' \
| awk 'NR==FNR{a[$1]=1;next}!($1 in a){print $0}' - "$parents"_marker \
> tmp && mv -f tmp "$parents"_marker

#the edge marker of contigs
awk 'NR==1{prev_chr=$1;first=$0;last=$0}NR>1{if($1!=prev_chr){print first;print last;first=$0}else{last=$0}prev_chr=$1}END{print first;print last}' "$parents"_marker > "$parents"_marker_edge
