module load GEMMA/0.98.1
module load BEDTools/2.27

mkdir gemma
mkdir SNPgemma
mkdir SVgemma
f=new_GT/ls361
kin=GT/gemma_Kinship/result.sXX.txt
q=GT/ls361_gemma_cov_pc3.txt
out=ls358_SNV_SV
p=$1
#gemma -bfile $f -k $kin -lmm 1 -c $q -o ${out}_${p} -p gemma_phes/${p}.txt
#grep "SNV" output/${out}_${p}.assoc.txt >SNPgemma/${out}_${p}.association.txt
#grep -v "SNV" output/${out}_${p}.assoc.txt |grep -v 'p_wald'>SVgemma/${out}_${p}.association.txt
#python gemma_sig_bed.py SVgemma/${out}_${p}.association.txt SVgemma/${out}_${p}_sig.bed $p
python gemma_sig_SNP_bed.py SNPgemma/${out}_${p}.association.txt SNPgemma/${out}_${p}_sig.bed $p
#bedtools closest -b /lguo/QYYang/zqyang/youcai_SV/ref/zs11.GLS.vGLS20211229.bed -a SVgemma/${out}_${p}_sig.bed -d|awk '{if($13<50000 && $13>=0) print $0}' |cut -f1-6,8-10,12,13 >SVgemma/${out}_${p}_sig_gene.txt
#python merge_gemma_dist_PVE.py SVgemma/${out}_${p}.association.txt SVgemma/${out}_${p}_sig_gene.txt $p SVgemma/${out}_${p}_sig_gene_PVE.txt
bedtools closest -b /data4/QYYang/zqyang/youcai_SV/ref/zs11.GLS.vGLS20211229.bed -a SNPgemma/${out}_${p}_sig.bed -d|awk '{if($13<50000 && $13>=0) print $0}' |cut -f1-6,8-10,12,13 >SVgemma/${out}_${p}_sig_gene.txt
python merge_gemma_dist_PVE.py SVgemma/${out}_${p}.association.txt SNPgemma/${out}_${p}_sig_gene.txt $p SNPgemma/${out}_${p}_sig_gene_PVE.txt

#python plot_mht_chr.py $p &
#python plot_mht.py $p
#python plot_fig_change_color.py -p emmax/${out}_${p}.ps -l ${out}_${p} -o emmax/${out}_${p}
#python plot_fig_change_color.py -p snp_emmax/gl505snp_${p}.ps -l 505SNP_${p} -o snp_emmax/gl505snp_${p} &
