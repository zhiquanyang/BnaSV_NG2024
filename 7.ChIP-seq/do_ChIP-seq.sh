module load fastp/0.23.0
module load Bowtie2/2.4.4
module load SAMtools/1.9
module load deepTools/3.5.0

ref=ref/zs11.genome.fa
i=$1
s=$2
fq1=cleandata/${i}_R1_001.fastq.gz
fq2=cleandata/${i}_R2_001.fastq.gz
cpu=4
if [ ! -d 01.trim ]
then
        mkdir 01.trim
fi

if [ ! -d 02.bam ]
then
    mkdir 02.bam
fi

if [ ! -d 03.macs ]
then
    mkdir 03.macs
fi

if [ ! -d 04.bw ]
then
    mkdir 04.bw
fi

#QC using fastp
fastp -i $fq1  -o  01.trim/${i}_1.fq.gz -I $fq2 -O 01.trim/${i}_2.fq.gz -j 01.trim/${i}_report.json -h 01.trim/${i}_json.html -w $cpu > 01.trim/${i}.trim.log

#Bowtie
bowtie2 -x $ref -p $cpu -1 $fq1 -2 $fq2 |samtools view -q 10 -bhS -o 02.bam/${i}.raw.bam

samtools flagstat 02.bam/${i}.raw.bam > 02.bam/${i}.bowtie2.log

samtools view -@ $cpu -bh -o 02.bam/${i}.bam -q 10 02.bam/${i}.raw.bam

samtools sort -@ $cpu -o 02.bam/${i}.sort.bam 02.bam/${i}.bam
samtools index 02.bam/${i}.sort.bam && rm 02.bam/${i}.raw.bam 02.bam/${i}.bam
java -jar ~/software/picard.jar MarkDuplicates \
I=02.bam/${i}.sort.bam \
O=02.bam/${i}.dd.bam \
M=02.bam/${i}_dup_metrics.txt \
REMOVE_DUPLICATES=true
samtools index 02.bam/${i}.dd.bam