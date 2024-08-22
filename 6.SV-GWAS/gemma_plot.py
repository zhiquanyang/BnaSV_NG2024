# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 10:27:24 2021

@author: zqyang

Organization: HuaZhong Agricutural University
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

infile=sys.argv[1]
outfile=sys.argv[2]
#infile=r"F:\cotton\cottonMD\GWAS\Du1782\Gemma\blup.FM.assoc.txt"
#outfile=r"F:\cotton\cottonMD\GWAS\Du1782\Gemma\blup.FM.assoc.png"
faifile="/home/zqyang/ref_genome/ZS11_new/zs11.genome.fa.fai"
title=sys.argv[3]
chrs=[]
pos=[]
ps=[]
chr2start={}
starts=[]
chrs=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','C01','C02',
      'C03','C04','C05','C06','C07','C08','C09']
nsnp=54885
cutoff=-np.log10(1/float(nsnp))
k=0
for line in open(faifile):
    lines=line.strip().split('\t')
    chri=lines[0].replace('scaffold','')
    leni=int(lines[1])+1000000
    if chri in chrs:
        if k==0:
            chr2start.setdefault(chri,0)
            s=leni
            starts.append(0)
            k+=1
        else:
            starts.append(s)
            chr2start.setdefault(chri,s)
            s=s+leni

plt.rcParams['font.sans-serif'] = 'Arial'      
fig=plt.figure(figsize=(8,4))
ax=plt.subplot(111)


#plot SNV
color_list=["#6B8E23","#FFA500"]
colors=[]
k=0
ps=[]

#plot SV
color_list=["#6F99AD","#FFDC91"]
colors=[]
k=0
ps=[]
for line in open(infile):
    lines=line.strip().split('\t')
    chri=lines[0]
    if chri=='chr':
        continue
    else:
        chri=int(chri)
        if(chri<=19):
            chri=chrs[chri-1]
    if chri in chrs:
        idi=lines[1]
        pos=int(lines[2])
        p=float(lines[-1])
        ps.append(p)
        if k==0:
            x=[]
            y=[]
            x.append(chr2start[chri]+pos)
            y.append(-np.log10(p))
            color=color_list[int(chrs.index(chri)%(len(color_list)))]
            k+=1
            chr0=chri
            ymax=np.max(np.array(y))*1.25
        elif chri==chr0:
            x.append(chr2start[chri]+pos)
            y.append(-np.log10(p))
        else:
            x=np.array(x)
            y=np.array(y)
            ymax=np.max([ymax,np.max(y*1.25)])
            x1=x[y<=cutoff]
            y1=y[y<=cutoff]
            ax.scatter(x1,y1,color=color,s=2,alpha=1)
            x2=x[y>cutoff]
            y2=y[y>cutoff]
            if x2.shape[0]>0:
                ax.scatter(x2,y2,color="red",s=2,alpha=1)
            x=[]
            y=[]
            x.append(chr2start[chri]+pos)
            y.append(-np.log10(p))
            color=color_list[int(chrs.index(chri)%(len(color_list)))]
            chr0=chri
x=np.array(x)
y=np.array(y)
ymax=np.max([ymax,np.max(y*1.25)])
x1=x[y<=cutoff]
y1=y[y<=cutoff]
ax.scatter(x1,y1,color=color,s=2,alpha=1)
x2=x[y>cutoff]
y2=y[y>cutoff]
if x2.shape[0]>0:
    ax.scatter(x2,y2,color="red",s=4,alpha=1)
ymax2=np.max(y*1.25)
xmin=0
xmax=s
ax.hlines(cutoff,xmin,xmax,linestyles="--", \
  label="cutoff={0:.2f}".format(cutoff),linewidth=1,colors='black')

plt.legend(edgecolor="none")
ax.set_xticks(starts)
ax.set_xticklabels(chrs,rotation=30)
ax.set_xlim(xmin,xmax)
ax.set_ylim(0,ymax)
ax.set_ylabel("-Log10(P-value)",fontsize="x-large")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig(outfile,dpi=600)
