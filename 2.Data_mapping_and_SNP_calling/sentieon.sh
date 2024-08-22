module load SAMtools/1.9
module load BWA/0.7.17
module load sentieon/201808.08
export SENTIEON_LICENSE=mn01:9000
fq1=$1
fq2=$2
i=$3
fasta=$4
nt=$5
mq=10
group_prefix='bna2105'
platform='illumina'
if [ ! -d 02.FinalBam ];then mkdir -p 02.FinalBam;fi
if [ ! -d 03.OriGVCF ];then mkdir -p 03.OriGVCF;fi
(bwa mem -M -R "@RG\tID:${group_prefix}_${i}\tSM:${i}\tPL:$platform" -t $nt -K 10000000 $fasta $fq1 $fq2 || echo -n 'error' ) |samtools view -o tmp/${i}_q0.bam -bh -
samtools flagstat tmp/${i}_q0.bam > 02.FinalBam/${i}.flagstat
samtools stats -@ $nt tmp/${i}_q0.bam > 02.FinalBam/${i}.stats
samtools view -h -q $mq tmp/${i}_q0.bam | sentieon util sort -r $fasta -o tmp/${i}_sorted.bam -t $nt --sam2bam -i -
samtools sort -@ $CPU -o tmp/${i}_sorted.bam tmp/${i}_q0.bam && rm tmp/${i}_q0.bam
java -jar ~/software/picard.jar MarkDuplicates \
I=tmp/${i}_sorted.bam \
O=02.FinalBam/${i}.bam \
M=02.FinalBam/${i}_dup_metrics.txt \
REMOVE_DUPLICATES=true

samtools idxstats 02.FinalBam/${i}.bam >02.FinalBam/${i}.idxstats
java -jar ~/software/picard.jar CollectWgsMetrics I=02.FinalBam/${i}.realn.bam O=tmp/${i}.CollectWgsMetrics R=${fasta}
sentieon driver -r $fasta -t $nt -i 02.FinalBam/${i}.bam --algo QualCal tmp/${i}_recal_data.table
sentieon driver -r $fasta -t $nt -i 02.FinalBam/${i}.bam -q tmp/${i}_recal_data.table --algo QualCal tmp/${i}_recal_data.table.post
sentieon driver -t $nt --algo QualCal --plot --before tmp/${i}_recal_data.table --after tmp/${i}_recal_data.table.post tmp/${i}_recal.csv
sentieon plot bqsr -o tmp/${i}_recal_plots.pdf tmp/${i}_recal.csv
sentieon driver -r $fasta -t $nt -i 02.FinalBam/${i}.bam -q tmp/${i}_recal_data.table --algo Haplotyper --emit_mode gvcf 03.OriGVCF/${i}_gvcf.gz