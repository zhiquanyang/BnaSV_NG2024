# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:43:12 2020

@author: zqyang

Organization: HuaZhong Agricultural University
"""

import numpy as np
import argparse

Args=argparse.ArgumentParser(description="SV filter step1")
Args.add_argument('-i','--input',dest='input')
Args.add_argument('-l','--len',dest='len',type=int,default=50)
Args.add_argument('-n','--nreads',dest='nreads',type=int,default=4)
Args.add_argument('-o','--out',dest='out')
args=Args.parse_args()

infile=args.input
outfile=args.out
cutoff=args.len
out=open(outfile,'w')
n=0
n_re=args.nreads
for line in open(infile):
    if line.startswith('#'):
       out.write(line)
    elif("UNRESOLVED" in line)or("IMPRECISE" in line):
        n+=1
        continue
    else:
        lines=line.strip().split('\t')
        try:
            r=True
            info_dict={}
            infos=lines[7].split(';')
            gt=lines[9].split(":")[0]
            if("1/1" not in gt):
                r=False
                n+=1
            else:
                for info_i in infos:
                    if("=" in info_i):
                        info_is=info_i.split('=')
                        info_dict.setdefault(info_is[0],info_is[1])
                
                if('[' in lines[4])or(']' in lines[4]):
                    sv_len=1
                else:
                    sv_len=np.abs(int(info_dict['SVLEN']))
                if sv_len<cutoff:
                    n+=1
                    r=False
                reads=np.abs(info_dict['RE'])
                if reads<n_re:
                    n+=1
                    r=False
            if r:
                out.write(line)
        except:
            out.write(line)
print("filter out {0} SVs".format(n))

