# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:47:57 2020

@author: zqyang

Organization: HuaZhong Agricutural University
"""

import sys
import numpy as np
infile=sys.argv[1]
outfile=sys.argv[2]

out=open(outfile,'w')
for line in open(infile):
    if(line.startswith("#")):
        out.write(line)
    else:
        lines=line.strip().split('\t')
        chri=lines[0]
        pos=int(lines[1])
        infos_i=lines[7]
        infos=infos_i.split(";")
        info_dict={}
        for info_i in infos:
            if("=" in info_i):
                info_is=info_i.split("=")
                info_dict.setdefault(info_is[0],info_is[1])
        ti=info_dict['SVTYPE']
        e=int(info_dict['END'])
        leni=np.abs(int(info_dict['SVLEN']))
        lines[2]="{0}_{1}_{2}_{3}_{4}".format(chri,pos,e,ti,leni)
        out.write("\t".join(lines)+'\n')
out.close()
