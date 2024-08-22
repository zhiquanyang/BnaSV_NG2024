# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 23:40:08 2019

@author: zqyang

Organization: HuaZhong Agricultural University
"""

import sys
infile=sys.argv[1]
outfile=sys.argv[2]
out=open(outfile,'w')
files=[line.strip() for line in open(infile)]
k=0
for infile in files:
    for line in open(infile):
        if line.startswith('#'):
           if k==0:
               out.write(line)
        else:
           out.write(line)
    k+=1
out.close()
