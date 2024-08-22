module load MACS2/2.2.7.1
module load ucsc_kentUtils/v389
treat=$1
control=$2
n=$treat

############MACS2 call peaks

g=800000000
#narrow peak for ChIP-seq
macs2 callpeak -t 02.bam/${treat}.bam -c 02.bam/${control}.bam -f BAM -n $n -B -g $g --outdir 03.macs/

macs2 bdgcmp -t 03.macs/${n}_treat_pileup.bdg -c 03.macs/${n}_control_lambda.bdg --outdir 03.macs/ -o ${n}_PE.bdg  -m FE
sort -k1,1 -k2,2n 03.macs/${n}_PE.bdg > 03.macs/${n}_PE.sort.bdg
bedGraphToBigWig 03.macs/${n}_PE.sort.bdg ref/zs11.genome.fa.fai 03.macs/${n}.bw && rm 03.macs/${n}_PE.bdg 03.macs/${n}_PE.sort.bdg
