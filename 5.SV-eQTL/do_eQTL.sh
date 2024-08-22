module load GSL/2.4
module load GEMMA/0.98.1
unset arr
declare -A arr
arr=( ['scaffoldA01']=1 ['scaffoldA02']=2 ['scaffoldA03']=3 ['scaffoldA04']=4 ['scaffoldA05']=5 ['scaffoldA06']=6 ['scaffoldA07']=7 ['scaffoldA08']=8 ['scaffoldA09']=9 ['scaffoldA10']=10 ['scaffoldC01']=11 ['scaffoldC02']=12 ['scaffoldC03']=13 ['scaffoldC04']=14 ['scaffoldC05']=15 ['scaffoldC06']=16 ['scaffoldC07']=17 ['scaffoldC08']=18 ['scaffoldC09']=19 )
infile=$1
f=GT/gl505_20DAF_SV
kin=GT/output/gl505_20DAF.sXX.txt
q=GT/gl505_20DAF_pca3_cov.txt
out=gl505_20DAF_SV
cat expressed_gene.bed|while read line
do
array=(${line})
chri0=${array[0]}
chri=${arr[$chri0]}
s=$((array[1]-1000000))
e=$((array[2]+1000000))
p=${array[3]}
phe_f=gene_exp/${p}.txt
gemma -bfile $f -k $kin -lmm 1 -p $phe_f -c $q -o ${p}
awk '{if($12<1/55634) print $0}' output/${p}.assoc.txt > output/${p}.sig_assoc.txt && awk '{if($1=="'$chri'"&&$3>="'$s'"+0&&$3<="'$e'"+0) print $0}' output/${p}.assoc.txt > output/${p}.cis_assoc.txt && rm output/${p}.assoc.txt
python gemma2pve.py output/${p}.sig_assoc.txt $p output/${p}.sig_pve_assoc.txt
done
