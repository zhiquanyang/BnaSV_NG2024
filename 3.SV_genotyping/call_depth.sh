s=$1
ref=/public/home/zqyang/youcai_SV/ref/zs11.genome.fa
bam=/public/home/zqyang/youcai_SV/04.pop/bams/${s}.bam
out=04.depth/01.depth
if [ ! -d $out ];then mkdir -p $out;fi
samtools depth -a --reference $ref $bam |awk '{ sum += $3; } END { print "'"$s"'" "\t" sum "\t" sum/NR }' >${out}/${s}.depth