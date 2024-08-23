# 5 SV-GWAS

Filtering SNPs, InDels and SVs before GWAS

```shell
vcf=/home/zqyang/data4/youcai_SV/04.pop/2105_maf0.01_miss0.7_SNP_SV_C02_GTR2_SV_edit.phase.vcf.gz
#vcf=/home/zqyang/data4/youcai_SV/04.pop/C02_GTR2/2105_maf0.01_miss0.7_SNP_SV.phase.vcf.gz
id_f=361inv.txt
out=bna361
vcftools --gzvcf $vcf --keep $id_f --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ${out}_maf0.05
python het_filter.py -i ${out}_maf0.05.recode.vcf -m 0.05 -c 0.5 -o ${out}_maf0.05.het0.5
zgrep -v "#" ${out}_maf0.05.het0.5.vcf.gz |cut -f3 |grep SNV > ${out}_maf0.05.het0.5_SNVids.txt
zgrep -v "#" ${out}_maf0.05.het0.5.vcf.gz |cut -f3 |grep -v SNV > ${out}_maf0.05.het0.5_SVids.txt
vcftools --gzvcf ${out}_maf0.05.het0.5.vcf.gz --snps ${out}_maf0.05.het0.5_SNVids.txt --recode --recode-INFO-all --out ${out/_SNV_SV/}_maf0.05_SNV
vcftools --gzvcf ${out}_maf0.05.het0.5.vcf.gz --snps ${out}_maf0.05.het0.5_SVids.txt --recode --recode-INFO-all --out ${out/_SNV_SV/}_maf0.05_SV
plink --vcf ${out}_maf0.05.het0.5.vcf.gz --make-bed --out ${out}

out=bna361
plink --vcf ${out}_maf0.05_SV.recode.vcf  --make-bed --out ${out}_SV
```

Generate kinship matrix and PCA results

```shell
gcta64 --bfile ls361 --autosome --make-grm --out ls361
gcta64 --grm ls361 --pca 10 --out ls361_pca10
awk '{print "1\t"$3"\t"$4"\t"$5}' ls361_pca10.eigenvec >ls361_pca3_cov.txt
cut -f2 ls361_pca3_cov.txt > tmp.txt
gemma -bfile ls361 -gk 2 -o ls361 -p tmp.txt
```

Perform GWAS using GEMMA.

```shell
#!/bin/bash
cpu=1
cat phes.txt|while read p
do
sh do_gemma_GWAS.sh $p
done
```
