# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 19:58:54 2019

@author: zqyang

Organization: HuaZhong Agricultural University
"""

import argparse
import pandas as pd
import numpy as np
Args=argparse.ArgumentParser(description="filter out FP SVs")
Args.add_argument('-i','--input',dest='input')
Args.add_argument('-b','--bed',dest='bed')
Args.add_argument('-d','--delta',dest='delta',type=int,help='|a-b|<=delta (10)bp',default=10)
Args.add_argument('-r','--ratio',dest='ratio',type=float,help='overlap ratio (0.8)',default=0.8)
Args.add_argument('-o','--out',dest='out')
args=Args.parse_args()

infile=args.input
bedfile=args.bed
out=open(args.out,'w')
d=args.delta
r=args.ratio

bed=pd.read_csv(bedfile,sep='\t',names=['chr','start','end','id'])
dat=pd.read_csv(infile,names=['chr1','s1','e1','id1','chr2','s2','e2','id2'],sep='\t')
allids=list(bed.loc[:,'id'])
filter_out=[]
for i in dat.index:
    s1,e1,id1,s2,e2=dat.loc[i,['s1','e1','id1','s2','e2']]
    len1=float(e1-s1+1)
    len2=float(e2-s2+1)
    ov=np.min([e1,e2])-np.max([s1,s2])+1
    ov1=ov/len1
    ov2=ov/len2
    if(np.abs(s1-s2)<=d)and(np.abs(e1-e2)<=d)and(ov1>=r)and(ov2>=r):
        filter_out.append(id1)
    else:
        continue
outids=sorted(list(set(allids)-set(filter_out)))

for i in outids:
    out.write(str(i)+"\n")
out.close()
    
