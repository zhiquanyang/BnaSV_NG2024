ref=$1
gtf=$2
fq1=$3
fq2=$4
cpu=$5
s=$6
mkdir 01.bams
mkdir 02.stringtie/${s}
hisat2 -p $cpu --summary-file 01.bams/${s}_summary.txt --dta -x $ref -1 $fq1 -2 $fq2 | samtools view -Sbho 01.bams/${s}.bam - 1>01.bams/${s}.log 2>01.bams/${s}.err
samtools sort -@ $cpu -o 01.bams/${s}.srt.bam  01.bams/${s}.bam
samtools index 01.bams/${s}.srt.bam && rm 01.bams/${s}.bam

mkdir -p bedgraph
bedtools genomecov -bg -split -ibam 01.bams/${s}.srt.bam -g /public/home/yhu/zqyang/youcai_SV/leaf_RNAseq/genome_size.txt > bedgraph/${s}.bedgraph
sort -k1,1 -k2,2n bedgraph/${s}.bedgraph > bedgraph/${s}.sort.bedgraph
mkdir -p bigwig
bedGraphToBigWig bedgraph/${s}.sort.bedgraph /public/home/yhu/zqyang/youcai_SV/leaf_RNAseq/genome_size.txt bigwig/${s}.bw
stringtie 01.bams/${s}.srt.bam -p $cpu -G $gtf -e -B -o 02.stringtie/${s}/${s}.gtf -A 02.stringtie/${s}/${s}.tab -l ${s}
