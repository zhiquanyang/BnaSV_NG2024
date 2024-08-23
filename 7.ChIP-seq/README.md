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
