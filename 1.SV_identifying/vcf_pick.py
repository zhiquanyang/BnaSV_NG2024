# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 18:01:13 2019

@author: zqyang

Organization: HuaZhong Agricultural University
"""

import argparse

Args=argparse.ArgumentParser(description='filtering vcf by SV-ids')
Args.add_argument('-i','--input',dest='input')
Args.add_argument('-v','--vcf',dest='vcf')
Args.add_argument('-o','--out',dest='out')
args=Args.parse_args()
vcffile=args.vcf
infile=args.input
outfile=args.out
out=open(outfile,'w')
#info=pd.read_table(info_file)
ids=sorted([int(i.strip()) for i in open(infile)])
k=0
for line in open(vcffile):
   if line.startswith('#'):
       out.write(line)
   else:
       lines=line.strip().split('\t')
       try:
           if lines[2]==str(ids[k]):
               out.write(line)
               k+=1
           else:
               continue
       except:
           out.write(line)
out.close()
           
