#input=03.OriGVCF/gangan_map_zs11.minimap2.vcf
#out=04.Filtervcf/gangan_map_zs11.minimap2
vcf=$1
outpre=$2
python sv2type.py ${vcf} ${outpre}
for i in DEL INV INS TRA DUP
do
input=${outpre}_${i}.vcf
out=${outpre}_${i}
python vcf_filter_step1.py -i $input -o ${out}_flt1.vcf
python vcf_filter_step2.py -i ${out}_flt1.vcf -o ${out}_flt2
python vcf_pick.py -i ${out}_flt2_filter_ids.txt -v ${out}_flt1.vcf -o ${out}_flt2.vcf
python vcf2bed.py -i ${out}_flt2.vcf -o ${out}_flt2.bed
bedtools intersect -a ${out}_flt2.bed -b 04.Filtervcf/zs11_zs11.ngmlr_${i}_flt2.bed -wa -wb >${out}_flt2_merge.txt
python vcf_filter3.py -i ${out}_flt2_merge.txt -o ${out}_flt3.txt -b ${out}_flt2.bed -d 1000000
python vcf_pick.py -i ${out}_flt3.txt -v ${out}_flt1.vcf -o ${out}_final.vcf
done
