# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 10:17:15 2019

@author: zqyang

Organization: HuaZhong Agricultural University
"""

import sys
infile=sys.argv[1]
outpre=sys.argv[2]
dels=open("{0}_DEL.vcf".format(outpre),'w')
invs=open("{0}_INV.vcf".format(outpre),'w')
ins=open("{0}_INS.vcf".format(outpre),'w')
tras=open("{0}_TRA.vcf".format(outpre),'w')
dups=open("{0}_DUP.vcf".format(outpre),'w')
#invdups=open("{0}_INVDUP.vcf".format(outpre),'w')
for line in open(infile):
   if line.startswith('#'):
       dels.write(line)
       invs.write(line)
       ins.write(line)
       tras.write(line)
       dups.write(line)
#       invdups.write(line)
   else:
       lines=line.strip().split('\t')
       infos=lines[7].split(';')
       info_dict={}
       for i in infos:
           info_i=i.split('=')
           if len(info_i)!=2:
               continue
           else:
               ii,ij=info_i[:2]
           if ii=='SVTYPE':
               if ij=="DEL":
                   dels.write(line)
               elif ij=="INV":
                   invs.write(line)
               elif ij=="INS":
                   ins.write(line)
               elif(ij=="BND")or(ij=="TRA"):
                   tras.write(line)
               elif ij=="DUP":
                   dups.write(line)
               elif ij=="INVDUP":
                   invs.write(line)
dels.close()
invs.close()
ins.close()
tras.close()
dups.close()
