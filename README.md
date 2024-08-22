# BnaSV_NG2024
# 1. SV identifying

## &#x20;1.1 Genome comparison

We compared the contigs between the ZS11 reference genome and the other 15 assembled genomes using NUCmer program in MUMmer4.

```shell
sh cpu=40
for i in zs11
do
        ref="genomes/${i}.genome.fa"
        for q in ganganF73 no2127 quintaA shengli3 tapidor3  westar zheyou73  zs11
        do
        qfile=contigs/${q}.contigs.fa
        pre=genome_SV_${q}_${i}
        time sh mummer_mum.sh $ref $qfile $i $q
        done
done

```

Merge SVs of all accessions to create the final one SV file of genome alignment

```shell
ls vcfs/*.vcf >nucdiff_vcf.txt
out=18samples_nucdiff.vcf
SURVIVOR merge nucdiff_vcf.txt 100 -1 1 1 -1 -1 $out
python add_seq_survivor.py nucdiff_vcf.txt $out 18samples_nucdiff_seq.vcf

```

## 1.2 Indentifying SVs by alignment of long reads

We used NGMLR (version 0.2.8)^14^ to map the long-read (> 500 bp) of each accession onto the ZS11 reference genome and carried out the SV calling by Sniffles (version 1.0.7)^14^.&#x20;

```shell
ref=zs11.genome.fa # set the reference genome
for s in gangan no2127 quinta shengli tapidor westar zheyou7 # set the long-read sequencing dataset of queried sample
do
q=rawdata/{s}_L500.fa
cpu=30
i={i}_ngmlr
seqtk seq -L 500 -a rawdata/{s}.fq.gz > rawdata/{s}_L500.fa # filter out reads of less than 500 bp
time sh ngmlr.sh $ref $q $i $i $cpu # 
done
```

## 1.3 Filtering SVs identified by long reads

We filtered the SVs with the MD flag. All SVs were filtered out with low-quality (flag: UNRESOLVED), ambiguous breakpoints (flag: IMPRECISE), less than four supporting reads, < 50 bp and duplicate calling.

```shell
sh vcf_filter.sh 03.OriGVCF/zs11_map_zs11.ngmlr.vcf 04.Filtervcf/zs11_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/tapidor_map_zs11.ngmlr.vcf 04.Filtervcf/tapidor_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/westar_map_zs11.ngmlr.vcf 04.Filtervcf/westar_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/gangan_map_zs11.ngmlr.vcf 04.Filtervcf/gangan_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/quinta_map_zs11.ngmlr.vcf 04.Filtervcf/quinta_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/shengli_map_zs11.ngmlr.vcf 04.Filtervcf/shengli_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/zheyou7_map_zs11.ngmlr.vcf 04.Filtervcf/zheyou7_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/no2127_map_zs11.ngmlr.vcf 04.Filtervcf/no2127_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/NY7_map_zs11.ngmlr.vcf 04.Filtervcf/NY7_zs11.ngmlr
sh vcf_filter.sh 03.OriGVCF/Express_617_map_zs11.ngmlr.vcf 04.Filtervcf/Express_617_map_zs11.ngmlr
```

Next, we merged SV sets from different accessions.

```shell
for i in DEL INS INV DUP TRA
do
echo $i
ls 04.Filtervcf/*_zs11.ngmlr_${i}_final.vcf > tmp/mapzs11_${i}_files_ngmlr.txt
SURVIVOR merge tmp/mapzs11_${i}_files_ngmlr.txt 10 1 1 1 0 50 04.Filtervcf/mapzs11_${i}_ngmlr.vcf
grep -v "#" 04.Filtervcf/mapzs11_${i}_ngmlr.vcf |wc -l
done
ls 04.Filtervcf/mapzs11_*_ngmlr.vcf >vcf_files_ngmlr.txt
python merge_vcf.py vcf_files_ngmlr.txt mapzs11_ngmlr.vcf
java -jar ~/software/picard.jar SortVcf I=mapzs11_ngmlr.vcf O=mapzs11_ngmlr.sort.vcf
```

Then, we merged SV sets from these two identification approaches using SURVIVOR (version 1.0.3) with the parameters “10 1 1 1 0 50”^68^, and combined them into a single VCF using the population calling method of the Sniffles pipeline^14^.

```shell
python do_vcf.py 18samples_nucdiff_seq.vcf contig.vcf contig
python do_vcf.py mapzs11_ngmlr.sort.vcf pacbio.vcf pacbio
ls contig.vcf pacbio.vcf > genome_pacbio_vcf_files.txt
SURVIVOR merge genome_pacbio_vcf_files.txt 100 1 1 -1 -1 50 genome_pacbio.vcf #delete some SVs just like SVTYPE="INS,DEL", "INS,INV"
python add_seq_merge_survivor.py genome_pacbio_vcf_files.txt genome_pacbio.vcf genome_pacbio2.vcf
python for_paragraph.py genome_pacbio2.vcf genome_pacbio_for_paragraph.vcf
```



# 2. Alignment of short reads from 2105 accessions and calling SNPs and InDels

## 2.1 Build the genome index

```shell
bwa index /public/home/zqyang/ref_genome/zs11.genome.fa
java -jar ~/software/picard.jar CreateSequenceDictionary R=/public/home/zqyang/ref_genome/zs11.genome.fa SP=ZS11 O=/public/home/zqyang/ref_genome/zs11.genome.dict
samtools faidx /public/home/zqyang/ref_genome/zs11.genome.fa
```

## 2.2 Mapping and calling SNPs and InDels

```shell
ref=/public/home/zqyang/ref_genome/zs11.genome.fa
cpu=30
cat samples.txt |while read i
cat ids.txt|while read i
do
bsub -J $i -n $cpu -R span[hosts=1] -o logs/${i}.out -e logs/${i}.err -q normal "sh sentieon.sh rawdata/${i}_1.fq.gz rawdata/${i}_2.fq.gz ${i} ${ref} ${cpu}"
done
```



# 3. SV genotyping for 2105 accessions

## 3.1 Calculate average sequencing depth for all accessions

```shell
#!/bin/bash
cat samples.txt | while read i;do
bsub -J $i -n 1 -o logs/${i}_depth.out -e logs/${i}_depth.err -q "normal" -R span[hosts=1] "sh call_depth.sh $i" done
```

## 3.2 SV genotyping using paragraph

```shell
#!/bin/bash
nt=4
ref=/public/home/zqyang/youcai_SV/ref/zs11.genome.fa
vcf=/public/home/zqyang/youcai_SV/02.smart/04.mergeSV/map_zs11_ngmlr_paragraph.vcf
if [ ! -d logs ];then mkdir -p logs;fi
#cat samples.txt | while read i;do
cat undo.txt|while read i;do
bsub -J ${i}_para -n $nt -o logs/${i}_paragraph.out -e logs/${i}_paragraph.err -q "normal" -R span[hosts=1] \
"
python run_paragraph.py -v $vcf -b bams -d 04.depth/01.depth/${i}.depth -r ${ref} -t $nt -o paragraph
"
done
```



# 4 RNA-seq

## 4.1 Create genome index

```shell
hisat2-build /data4/youcai_SV/ref/zs11.genome.fa /data4/youcai_SV/ref/zs11
```

## 4.2 Mapping and calculate gene expression levels

```shell
ref=/data4/youcai_SV/ref/zs11
fa=/data4/youcai_SV/ref/zs11.genome.fa
cpu=2
bsub -J ${sample} -n $cpu -q normal -o logs/${sample}.out -e logs/${sample}.err -R span[hosts=1] "hisat2_stringtie.sh $ref /data4/youcai_SV/ref/zs11.v0.gtf rawdata/${sample}_1.clean.fq.gz rawdata/${sample}_2.clean.fq.gz 9 ${sample}"
done
```

## 4.3 Merge gene expression levels of all accessions

```shell
python merge_tpm.py samples.txt all_sample_tpm.txt
```

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

# 6 SV-eQTL

For all expressed genes, perform association analysis using SV genotypes and gene expression levels in population

```shell
sh do_eQTL.sh
```

# 7. ChIP-seq data analysis

## 7.1 Build the genome index for reference genome

```shell
bowtie2-build --threads 10 /data4/youcai_SV/ref/zs11.genome.fa /data4/youcai_SV/ref/zs11
```

7.2 Mapping all reads of sequenced samples

First, prepare a file about information of all alibraries, just like:

```bash
IN_795-2_S19_L002_filtered      ZS11-1S_IN
INPUT_795-1_S12_L002_filtered   ZS11-1P_IN
INPUT_795-3_S13_L002_filtered   ZS11-2P_IN
INPUT_795-4_S14_L002_filtered   ZS11-2S_IN
INPUT_795-5_S15_L002_filtered   ZY821-1P_IN
INPUT_795-6_S16_L002_filtered   ZY821-1S_IN
INPUT_795-7_S17_L002_filtered   ZY821-2P_IN
INPUT_795-8_S18_L002_filtered   ZY821-2S_IN
IP_795-1_H3K27Ac_S5_L002_filtered       ZS11-1P_IP
IP_795-2_H3K27AC_S20_L002_filtered      ZS11-1S_IP
IP_795-3_H3K27Ac_S6_L002_filtered       ZS11-2P_IP
IP_795-4_H3K27Ac_S7_L002_filtered       ZS11-2S_IP
IP_795-5_H3K27Ac_S8_L002_filtered       ZY821-1P_IP
IP_795-6_H3K27Ac_S9_L002_filtered       ZY821-1S_IP
IP_795-7_H3K27Ac_S10_L002_filtered      ZY821-2P_IP
IP_795-8_H3K27Ac_S11_L002_filtered      ZY821-2S_IP
```

Next, map all reads of sequenced samples to the reference genome

```shell
cat sample_info.txt |while read i
do
        info=($i)
        id_i=${info[0]}
        s=${info[1]}
        bsub -J ${id_i} -n $cpu -R span[hosts=1] -o logs/${id_i}.out -e logs/${id_i}.err -q normal "time sh do_ChIP-seq.sh ${id_i} $s"
done
```

Finally, call peaks for case-control groups

```shell
cpu=1
cat compare_ids.txt|while read line
do
array=($line)
s=${array[0]}
i=${array[1]}
bsub -J ${s} -n $cpu -R span[hosts=1] -o logs/${s}.out -e logs/${s}.err -q normal " time sh callpeak.sh $s $i"
done
```

The compare\_ids.txt file is just like:

```bash
ZS11-1P_IP      ZS11-1P_IN
ZS11-2P_IP      ZS11-2P_IN
ZS11-1S_IP      ZS11-1S_IN
ZS11-2S_IP      ZS11-2S_IN
ZY821-1P_IP     ZY821-1P_IN
ZY821-2P_IP     ZY821-2P_IN
ZY821-1S_IP     ZY821-1S_IN
ZY821-2S_IP     ZY821-2S_IN
```

