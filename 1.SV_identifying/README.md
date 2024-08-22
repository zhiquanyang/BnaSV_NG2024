# 1. SV identifying

## &#x20;1.1 Genome alignment

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

