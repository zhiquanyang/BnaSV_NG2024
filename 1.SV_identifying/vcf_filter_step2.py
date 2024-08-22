# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:29:05 2019

@author: zqyang

Organization: HuaZhong Agricultural University


scaffoldA01     243384  9       TA       N       .       PASS    IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=scaffoldA01;END=243698;STD_quant_start=21.257940;STD_quant_stop=15.671630;Kurtosis_quant_start=-1.358364;Kurtosis_quant_stop=1.382319;SVTYPE=DEL;SUPTYPE=AL,SR;SVLEN=-314;STRANDS=+-;RE=12      GT:DR:DV        ./.:.:12
"""

import argparse
import numpy as np
import pandas as pd

Args=argparse.ArgumentParser(description="SV filter step1")
Args.add_argument('-i','--input',dest='input')
Args.add_argument('-o','--out',dest='out',help='we would create *_cluster_ids.txt, *_uniq_ids.txt and *_info.txt')
args=Args.parse_args()


infile=args.input
outpre=args.out

out1=open("{0}_cluster_ids.txt".format(outpre),'w')
out2=open("{0}_uniq_ids.txt".format(outpre),'w')
out3="{0}_info.txt".format(outpre)

print("step1: clustering SVs by positions ...")   
k=0
chrs=[]
starts=[]
ends=[]
ids=[]
ns=[]
for line in open(infile):
   if line.startswith('#'):
       continue
   else:
        lines=line.strip().split('\t')
        chri=lines[0]
        si=int(lines[1])
        idi=lines[2]
        if(('[' in lines[4])or(']' in lines[4]))and('SVTYPE=INV' not in lines[7]):
            ei=si+1
        for info in lines[7].split(';'):
            if(info.startswith('END='))and(not((('[' in lines[4])or(']' in lines[4]))and('SVTYPE=INV' not in lines[7]))):
                ei=int(info.split('=')[-1])
            if 'RE=' in info:
                ni=np.abs(int(info.split('=')[-1]))
        chrs.append(chri)
        starts.append(np.min([si,ei]))
        ends.append(np.max([si,ei]))
        ids.append(idi)
        ns.append(ni)
dat=pd.DataFrame({'chr':chrs,'start':starts,'end':ends,'id':ids,'n_reads':ns})
dat=dat.sort_values(by=['chr','start'])
dat.to_csv(out3,index=False,sep='\t')
ids=[]
#        out3.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chri,si,ei,idi,ni))
for i in dat.index:
    chri,si,ei,idi=dat.loc[i,['chr','start','end','id']]
    if k==0:
        ids=[idi]
        k+=1
        chr0=chri
        s0=si
        e0=ei
    elif chr0==chri:
        if(si>e0):
            if len(ids)>1:
                out1.write(",".join(ids)+"\n")
            else:
                out2.write(ids[0]+"\n")
            ids=[idi]
            s0=si
            e0=ei
            chr0=chri
        else:
            ids.append(idi)
            s0=np.min([s0,si])
            e0=np.max([e0,ei])
    else:
        if len(ids)>1:
            out1.write(",".join(ids)+"\n")
        else:
            out2.write(ids[0]+"\n")
        ids=[idi]
        s0=si
        e0=ei
        chr0=chri
if len(ids)>1:
    out1.write(",".join(ids)+"\n")
else:
    out2.write(ids[0]+"\n")
out1.close()
out2.close()
           
print("step2: pick uniq ids by supporting read number ...")
infile1="{0}_info.txt".format(outpre)
#info
infile2="{0}_cluster_ids.txt".format(outpre)
#cluter ids
infile3="{0}_uniq_ids.txt".format(outpre)
#uniq ids
outfile="{0}_filter_ids.txt".format(outpre)
out=open(outfile,'w')
#infile1=r"F:\油菜SV\02.SMART\gangan_map_zs11_flt2.minimap2_info.txt"
#infile2=r"F:\油菜SV\02.SMART\gangan_map_zs11_flt2.minimap2_cluster_ids.txt"
dat=pd.read_csv(infile1,sep='\t')
dat.index=list(dat.loc[:,'id'])
id_groups=[]
for line in open(infile2):
    id_groups.append([int(i) for i in line.strip().split(',')])

def filter_ids(dat,gi):
    ids_pass=[]
    dati=dat.loc[gi,]
    dati=dati.sort_values(by='n_reads',ascending=False)
    
    
    now_ids=list(dati.loc[:,'id'])
    tmp=now_ids[1:]
    
    id0=now_ids[0]
    ids_pass.append(id0)
    ids_left=[]
    s0,e0=dati.loc[id0,['start','end']]
    for i in tmp:
        si,ei=dati.loc[i,['start','end']]
        if(ei<s0)or(si>e0):
            ids_left.append(i)
    #初次筛选，取出最高reads支持，去掉重叠区域
    
    #correct wrong
    while len(ids_left)>1:
        my_ids_left=[]
        dat2=dat.loc[ids_left,]
        dat2=dat2.sort_values(by='start')
        k=0
        for my_i in dat2.index:
            my_s,my_e,myid=dat2.loc[my_i,['start','end','id']]
            if k==0:
                myids=[myid]
                k+=1
                my_s0=my_s
                my_e0=my_e
            else:
                if(my_s>my_e0):
                    if len(myids)>1:
                        my_ids_left=my_ids_left+myids
                    else:
                        ids_pass.append(myids[0])
                    myids=[myid]
                    my_s0=my_s
                    my_e0=my_e
                else:
                    myids.append(myid)
                    my_s0=np.min([my_s0,my_s])
                    my_e0=np.max([my_e0,my_e])
        if len(myids)>1:
            my_ids_left=my_ids_left+myids
        
        if len(my_ids_left)>1:
            dati=dat.loc[my_ids_left,]
            dati=dati.sort_values(by='n_reads',ascending=False)
            now_ids=list(dati.loc[:,'id'])
            tmp=now_ids[1:]
            id0=now_ids[0]
            ids_pass.append(id0)
            ids_left=[]
            s0,e0=dati.loc[id0,['start','end']]
            for i in tmp:
                si,ei=dati.loc[i,['start','end']]
                if(ei<s0)or(si>e0):
                    ids_left.append(i)
#            print(ids_left)
#            print("{0} left".format(len(ids_left)))
        else:
            ids_left=my_ids_left
    if len(ids_left)==1:
        ids_pass.append(ids_left[0])
    dati=dat.loc[ids_pass,:]
    dati=dati.sort_values(by='start')
    myids=list(dati.loc[:,'id'])
    return(myids)
out_ids=[]
for i in range(len(id_groups)):
    gi=id_groups[i]
    out_ids=out_ids+filter_ids(dat,gi)
uniqids=[int(line.strip()) for line in open(infile3)]
outids=sorted(out_ids+uniqids)
print("{0} left".format(len(outids)))
print("step3: output ids...")
for i in outids:
    out.write(str(i)+'\n')
out.close()
