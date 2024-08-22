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
