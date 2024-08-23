module load GEMMA/0.98.1
module load BEDTools/2.27

mkdir gemma
mkdir SNPgemma
mkdir SVgemma
f=new_GT/ls361
kin=GT/gemma_Kinship/result.sXX.txt
q=GT/ls361_gemma_cov_pc3.txt
out=ls358_SNV_SV
p=$1
gemma -bfile $f -k $kin -lmm 1 -c $q -o ${out}_${p} -p gemma_phes/${p}.txt
grep "SNV" output/${out}_${p}.assoc.txt >SNPgemma/${out}_${p}.association.txt
grep -v "SNV" output/${out}_${p}.assoc.txt |grep -v 'p_wald'>SVgemma/${out}_${p}.association.txt
