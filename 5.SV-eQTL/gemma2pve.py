# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 10:41:52 2021

@author: zqyang

Organization: HuaZhong Agricutural University
"""

import sys
import pandas as pd
import numpy as np

infile1=sys.argv[1]
outfile=sys.argv[3]
phe=sys.argv[2]

phe_f="/home/zqyang/data4/youcai_SV/12.Pathway/SV_eQTL_20220101/gl505_seed_20DAF/gl505_20DAF_seed_exp/{0}.txt".format(phe)
dat=pd.read_csv(infile1,sep='\t',names=['chr','SV','ps','n_miss','allele1',
                                        'allele0','af','beta','se','logl_H1',
                                        'l_remle','p_wald'])

values=[]
for line in open(phe_f):
    v=line.strip()
    if(v != 'NA'):
        values.append(float(v))
var_y=np.std(values)**2
print(str(var_y))
n=int(dat.shape[0])
if(n==0):
    out=open(outfile,'w')
    headers=['gene','chr','SV','ps','n_miss','allele1','allele0','af','beta','se',
               'logl_H1','l_remle','p_wald','r2']
    out.write('\t'.join(headers)+'\n')
    out.close()
else:
    #PVE=beta^2*2*maf*(1-maf)/var(y)
    dat.loc[:,'r2']=dat['beta']**2*2*dat['af']*(1-dat['af'])/var_y
    dat.loc[:,'gene']=phe
    dat=dat.loc[:,['gene','chr','SV','ps','n_miss','allele1','allele0','af','beta','se',
               'logl_H1','l_remle','p_wald','r2']]
    dat.to_csv(outfile,sep='\t',index=False)
