# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 16:11:28 2019

@author: zqyang

Organization: HuaZhong Agricultural University
"""

import argparse,os

Args=argparse.ArgumentParser(description="Run paragraph")
Args.add_argument('-v','--vcf',dest='vcf')
#Args.add_argument('-s','--sample',dest='sample')
Args.add_argument('-b','--bam',dest='bam')
Args.add_argument('-d','--depth',dest='depth')
Args.add_argument('-r','--ref',dest='ref')
Args.add_argument('-t','--cpu',dest='cpu',type=int)
Args.add_argument('-o','--out',dest='out')
args=Args.parse_args()
#multigrmpy.py -i my_ues.vcf -m 1000.txt -r ../ref/zs11.genome.fa -o 1000 --threads 5 -M 200 1>out.log 2>out.err &
if not os.path.exists(args.out):
    os.mkdir(args.out)
mdir=os.path.join(args.out,"MANIFEST")

if not os.path.exists(mdir):
    os.mkdir(mdir)
k=0
for line in open(args.depth):
    lines=line.strip().split('\t')
    if k==0:
        s=lines[0]
        dep=int(float(lines[-1])+0.5)
mfile=os.path.join(mdir,"{0}.txt".format(s))
bamfile=os.path.join(args.bam,"{0}.bam".format(s))
mout=open(mfile,'w')
mout.write("id\tpath\tdepth\tread length\n")
mout.write("{0}\t{1}\t{2}\t{3}\n".format(s,bamfile,dep,150))
mout.close()

outfile=os.path.join(args.out,s)
cmd="multigrmpy.py -i {0} -m {1} -r {2} -o {3} --threads {4} -M {5}".format(
        args.vcf,mfile,args.ref,outfile,args.cpu,dep*20)
print(cmd)
os.system(cmd)