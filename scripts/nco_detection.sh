#!/bin/sh

wd="path/to/wd" && cd $wd
vcf="path/to/vcf_files"               #GT of vcf files from offspring short reads VS parental genome

parents="P"           # P or M

awk 'NR==FNR{a[$1"_"$2]=1;next}$1"_"$2 in a{print $0}' "$parents"_marker $vcf \
| awk 'NR==FNR{a[$2"_"$3"_"$4"_"$5]=1;b[$7"_"$6"_"$5"_"$4]=1;next}(($1"_"$2"_"$4"_"$5 in a) || ($1"_"$2"_"$4"_"$5 in b)){print $0}' pvss_diff_"$parents" - \
| sed -e 's/\.\/././g' -e 's/0\/1/M/g' -e 's/0|1/M/g' -e 's/1\/1/1/g' -e 's/1|1/1/g' > "$parents"_conversion_var.txt

# NCO detection
awk '{print $1";"$2,$2}' "$parents"_marker > temp1_"$parents"
for i in $(seq 10 $(awk '$1!~"#"{print NF}' $vcf | sort | uniq));do
# romove contigs with SNP at the beginning and end
awk '{print $1,$2,$'$i'}' "$parents"_conversion_var.txt | awk '$3==1' \
| awk 'NR==FNR{a[$1]=1;next}!($1 in a){print $1";"$2}' <(awk '{print $1,$2,$'$i'}' "$parents"_conversion_var.txt | awk '$3==1' \
| awk 'NR==FNR{a[$1"_"$2]=1;next}$1"_"$2 in a{print $1}' "$parents"_marker_edge - | sort | uniq) - > temp2_"$parents"_"$i"
awk 'NR==FNR{a[$1]=1;next}$2";"$3 in a{print $2,$3,$4,$5}$7";"$6 in a{print $7,$6,$5,$4}' temp2_"$parents"_"$i" pvss_diff_"$parents" > "$parents"_"$i"_variants.txt
# merge adjacent markers
python ${wd}/conversion_detection.py temp1_"$parents" temp2_"$parents"_"$i" "$parents"_"$i"_conversion.txt
awk -F '[; ]' '{print $2,$1}' "$parents"_"$i"_conversion.txt | awk -F, '{print NF,$0}' > tmp && mv -f tmp "$parents"_"$i"_conversion.txt
# remove uncertain inherited contigs
# the file("$parents"_"$i"_inherit) contains inherited contig list
awk 'NR==FNR{a[$1]=1;next}$2 in a{print $0}' "$parents"_"$i"_inherit "$parents"_"$i"_conversion.txt > "$parents"_"$i"_conversion_filter.txt
awk 'NR==FNR{a[$1]=1;next}$2 in a{print $0}' "$parents"_"$i"_inherit "$parents"_"$i"_variants.txt > "$parents"_"$i"_variants_filter.txt
done && rm temp*

# remove 4 haplotypes
cat \
<(cat "$parents"_*_variants_filter.txt | sort | uniq | awk '$1~"primary"' | awk 'NR==FNR{a[$1"_"$2]=1;next}$2"_"$3 in a{print $7"_"$6}' - pvss_diff_"$parents") \
<(cat "$parents"_*_variants_filter.txt | sort | uniq | awk '$1~"secondary"' | awk 'NR==FNR{a[$1"_"$2]=1;next}$7"_"$6 in a{print $2"_"$3}' - pvss_diff_"$parents") \
> "$parents"_potential_co.txt
awk '{print $1";"$2,$2}' "$parents"_marker > temp1_"$parents"
for i in $(seq 10 $(awk '$1!~"#"{print NF}' $vcf | sort | uniq));do
awk 'NR==FNR{a[$1]=1;next}!($1"_"$2 in a){print $0}' "$parents"_potential_co.txt "$parents"_"$i"_variants_filter.txt > "$parents"_"$i"_variants_filter2.txt
awk '{print $1";"$2}' "$parents"_"$i"_variants_filter2.txt > temp2_"$parents"_"$i"
python ${wd}/conversion_detection.py temp1_"$parents" temp2_"$parents"_"$i" "$parents"_"$i"_conversion_filter2.txt
awk -F '[; ]' '{print $2,$1}' "$parents"_"$i"_conversion_filter2.txt | awk -F, '{print NF,$0}' > tmp && mv -f tmp "$parents"_"$i"_conversion_filter2.txt
done && rm temp*
