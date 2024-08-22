# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 13:49:29 2019

@author: zqyang

Organization: HuaZhong Agricultural University

scaffoldA01     NucDiff_v2.0    SO:0001873      501     501     .       .       .       ID=SV_114064.1;Name=translocation-insertion;ins_len=951;query_sequence=tapidor3v0_A03:39429705-39993194;query_coord=212445;breakpoint_query=212446-213396;blk_query=191229-212445;blk_ref=501-21719;blk_query_len=21217;blk_ref_len=21219;color=#A0A0A0
"""

import sys
import numpy as np
infile=sys.argv[1]
outfile=sys.argv[2]
s=sys.argv[3]
out=open(outfile,'w')
chrs=["scaffoldA01","scaffoldA02","scaffoldA03","scaffoldA04","scaffoldA05",
      "scaffoldA06","scaffoldA07","scaffoldA08","scaffoldA09","scaffoldA10",
      "scaffoldC01","scaffoldC02","scaffoldC03","scaffoldC04","scaffoldC05",
      "scaffoldC06","scaffoldC07","scaffoldC08","scaffoldC09"]

out=open(outfile,'w')
k=0
n=0
for line in open(infile):
   if line.startswith('##'):
       out.write(line)
   elif line.startswith('#'):
       lines=line.strip().split('\t')
       outlines=lines[:8]+['FORMAT',s]
       out.write('\t'.join(outlines)+'\n')
   else:
       lines=line.strip().split('\t')
       info_dict={}
       infoi=lines[7]
       infos=infoi.split(';')
       chri=lines[0]
       pos=int(lines[1])
       if chri not in chrs:
           continue
       info_ids=[]
       for i in infos:
           info_is=i.split('=')
           info_dict.setdefault(info_is[0],info_is[1])
           info_ids.append(info_is[0])
       ti=info_dict['SVTYPE']
       lines[2]="{0}_{1}_{2}_{3}_{4}".format(chri,pos,info_dict['END'],ti,np.abs(int(info_dict['SVLEN'])))
       outlines=lines[:8]+['GT','1/1']
       out.write('\t'.join(outlines)+'\n')
out.close()
