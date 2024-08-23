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
