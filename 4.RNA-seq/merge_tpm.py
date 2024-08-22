# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:25:11 2020

@author: zqyang

Organization: HuaZhong Agricultural University

"""

import sys,os,re
import pandas as pd

infile=sys.argv[1]
outfile=sys.argv[2]
samples=[i.strip() for i in open(infile)]
k=0

def read_tab(infile):
    genes=[]
    fpkm=[]
    for line in open(infile):
        if line.startswith('#'):
            continue
        else:
            lines=line.strip().split('\t')
            #print(line)
            if(line.startswith("Gene")):
                continue
            genes.append(lines[0])
            fpkm.append(float(lines[-1]))
    dat=pd.DataFrame({'Gene':genes,'TPM':fpkm})
    dat=dat.loc[:,['Gene','TPM']]
    return(dat)
for i in samples:
    infile=os.path.join("02.stringtie",i,"{0}.tab".format(i))
    dat=read_tab(infile)
    dat=dat.loc[:,['Gene',"TPM"]]
    dat.columns=['Gene',i]
    if(k==0):
        res=dat
        k+=1
    else:
        res=pd.merge(left=res,right=dat,on="Gene")
    print(i)
res.to_csv(outfile,sep='\t',index=False)
