# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 20:33:32 2020

@author: zqyang

Organization: HuaZhong Agricutural University
"""

import sys
import pandas as pd
infile=sys.argv[1]
merge_file=sys.argv[2]
outfile=sys.argv[3]
k=0
files=[line.strip() for line in open(infile)]
for filei in files:
    print(filei)
    ids=[]
    chrs=[]
    starts=[]
    ends=[]
    lens=[]
    seqs=[]
    myinfos=[]
    for line in open(filei):
        if(line.startswith("#")):
            continue
        else:
            lines=line.strip().split("\t")
            chri=lines[0]
            s=int(lines[1])
            infos_i=lines[7]
            infos=infos_i.split(";")
            info_dict={}
            for info_i in infos:
                if("=" in info_i):
                    info_is=info_i.split("=")
                    info_dict.setdefault(info_is[0],info_is[1])
            ti=info_dict['SVTYPE']
            e=int(info_dict['END'])
            leni=info_dict['SVLEN']
            #print(line)
            if(ti == "INS")or(ti == "DUP"):
                chrs.append(chri)
                starts.append(s)
                ends.append(e)
                lens.append(leni)
                ids.append(lines[2])
                myinfos.append(infos_i)
                seqs.append(info_dict['SEQ'])
    if(len(ids)==0):
        continue
    dati=pd.DataFrame({'id':ids,'chr':chrs,'start':starts,'end':ends,'info':myinfos,'seq':seqs,'len':lens})
    dati.loc[:,'source']=filei.replace("vcfs/","").replace("_zs11.vcf","")
    if(k==0):
        dat=dati
        k+=1
    else:
        dat=pd.concat([dat,dati])
dat.to_csv("ins_seq.txt",index=False,sep='\t')
out=open(outfile,'w')
drop=0
for line in open(merge_file):
    if(line.startswith('##')):
        out.write(line)
    elif(line.startswith('#')):
        lines=line.strip().split('\t')
        inds=lines[9:]
        out.write(line)
    else:
        lines=line.strip().split("\t")
        infos=lines[7].split(";")
        info_dict={}
        for info_i in infos:
            if("=" in info_i):
                info_is=info_i.split("=")
                info_dict.setdefault(info_is[0],info_is[1])
        ti=info_dict['SVTYPE']
        if(ti != "INS"):
            out.write(line)
        else:
            idi=lines[2]
            info_names=lines[8].split(':')
            k=info_names.index('ID')
            labels=[]
            for gti in lines[9:]:
                gtis=gti.split(':')
                labels.append(gtis[k])
            id_k=labels.index(idi)
            indi=inds[id_k]
            print(idi)
            dati=dat.loc[(dat['source']==indi)&(dat['id']==idi),:]
            s=list(dati['start'])[0]
            info_i=list(dati['info'])[0].replace("SVTYPE=DUP;","SVTYPE=INS;")
            lines[7]=info_i
            lines[1]=str(s)
            out.write('\t'.join(lines)+'\n')
out.close()
          
print("drop {0} seqs".format(drop))
