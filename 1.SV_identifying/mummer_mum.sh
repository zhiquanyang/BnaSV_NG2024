ref=$1
q=$2
rid=$3
qid=$4
#ref=../genomes/zs11.genome.fa
#q=contigs/wester.contigs.fa
cpu=50
pre="${qid}_${rid}_mum"
nucmer --mum --noextend -t $cpu -L 1000 -p $pre $ref $q
mkdir ${pre}
delta-filter -1 -i 95 ${pre}.delta >${pre}/${pre}.delta
python ~/software/NucDiff/nucdiff.py --proc $cpu --vcf yes --delta_file ${pre}/${pre}.delta $ref $q $pre $pre
