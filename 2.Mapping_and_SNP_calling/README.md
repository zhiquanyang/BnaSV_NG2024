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
