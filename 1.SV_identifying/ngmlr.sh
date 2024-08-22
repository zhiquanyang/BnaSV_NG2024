ref=$1
q=$2
out=$3
i=$4
cpu=20
pl="ONT"
#ngmlr -x ont -t $cpu -r $ref -q $q --rg-id $i --rg-sm $i --rg-lb $i --rg-pl $pl |samtools view -bhS -o ${out}.ngmlr.bam -
#sambamba view -t $cpu -o ${out}.ngmlr.bam ${out}.ngmlr.sam && rm ${out}.ngmlr.sam
#sambamba sort -o ${out}.ngmlr.sort.bam -t $cpu ${out}.ngmlr.bam && rm ${out}.ngmlr.bam
sniffles -t $cpu -m ${out}.ngmlr.sort.bam -v ${out}.ngmlr.vcf
