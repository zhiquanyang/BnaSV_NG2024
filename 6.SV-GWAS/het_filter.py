# -*- coding: utf-8 -*-
"""
Created on Tue Jan 02 16:49:43 2018

@author: zqyang

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	10-2	100-1
"""

import gzip
import argparse
import numpy as np

description='Filter vcf file according to Hybrid rate'
Args=argparse.ArgumentParser(description=description)
Args.add_argument('-i','--input',dest='input',type=str,help='input file')
Args.add_argument('-c','--het',dest='het',type=float,default=0.5)
Args.add_argument('-m','--maf',dest='maf',type=float,default=0.05)
Args.add_argument('-o','--prefix',dest='prefix',type=str,help='output file')
args=Args.parse_args()
infile=args.input
cutoff=args.het
maf=args.maf
if infile.endswith('gz'):
     filein=gzip.open(infile)
else:
     filein=open(infile)
#infile=r'F:\vcf_deal\data\test.vcf'
#out1=open(r'{0}.vcf'.format(sys.argv[2]),'w')
#out2=open(r'{0}.stat'.format(sys.argv[2]),'w')

#infile=r'F:\vcf_deal\data\test.vcf'
chrs=['scaffoldA01','scaffoldA02','scaffoldA03','scaffoldA04','scaffoldA05','scaffoldA06','scaffoldA07','scaffoldA08','scaffoldA09','scaffoldA10',
'scaffoldC01','scaffoldC02','scaffoldC03','scaffoldC04','scaffoldC05','scaffoldC06','scaffoldC07','scaffoldC08','scaffoldC09']
out1=gzip.open(r'{0}.vcf.gz'.format(args.prefix),'w')
out2=gzip.open(r'{0}.stat.gz'.format(args.prefix),'w')
out2.write("CHROM\tPOS\tID\tREF\tALT\tHET\tMAF\n")
gt_dict={'0/0':0,'0|0':0,'1/1':2,'1|1':2,'0/1':1,'0|1':1,'1/0':1,'1|0':1,\
'./.':-9,'.|.':-9}
for line in filein:
    if line[0] == '#':
        if line[0]=='##':
            continue
        else:
            out1.write(line)
            continue
    lines=line.strip().split('\t')
    chri=lines[0]
#    if(chri in chrs):
#        lines[0]=int(chrs.index(chri)+1)
#    else:
#        continue
    out_lines=[]
#    lines[2]="{0}_{1}_SNV".format(lines[0],lines[1])
    info=lines[:9]
    info[8]='GT'
    gts=[i.split(':')[0] for i in lines[9:]]
    out_lines=info+gts
    gt_num=[gt_dict[i] for i in gts]
    n=float(len(gt_num))
    het=gt_num.count(1)/n
    freq1=gt_num.count(0)/n
    freq2=gt_num.count(2)/n
    min_freq=np.min([freq1,freq2])
    out2_lines=info[:5]+[het,min_freq]
    out2.write('\t'.join(map(str,out2_lines))+'\n')
    if(het<=cutoff)and(min_freq>=maf)and(het<min_freq):
         out1.write('\t'.join(map(str,out_lines))+'\n')
out1.close()
out2.close()

