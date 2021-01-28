# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
### mapping should start from the first base pair of reference, which means the POS should be 1 of all reads
import os
from os import path
import sys
from Bio import pairwise2
from Bio.Seq import Seq 
import difflib


# %%
# indicate input files
#ins_sum='../../data/HQ4641_A_12D.res.txt'
#ref_file='../../data/VB18_F3_short.fa'
#out_put='/Users/lengsiewyeap/Desktop/code/Microhomology/results/insertion/HQ4641_A_12D.mh_res.txt'

ins_sum=sys.argv[1]
ref_file=sys.argv[2]
out_put=sys.argv[3]

# %%
def similar_equal(seq1,seq2):
    if len(seq1)==len(seq2):
        sim_len=0
        for i in range(len(seq1)):
            if seq1.lower()[i]==seq2.lower()[i]:
                sim_len+=1
        return(sim_len/len(seq1))
# %%
if not path.exists(out_put):
    os.makedirs(out_put)

# reference sequnce
ref_seq=''
for line2 in open(ref_file):
    if line2.strip() != "" and not line2.strip().startswith('>'):
        ref_seq=line2.strip()


# %%
f=open(ins_sum,'r')
h1=f.readline()
all_ins_microhomo=[]
for line in f:
    if int(line.strip().split('\t')[28]) in [1,2,3,11] and int(line.strip().split('\t')[20])>=20: # start position should bigger than 20
        tmp=line.strip().split('\t')
        ins_len=int(tmp[24])
        ins_seq=tmp[25]
        seq_start=int(tmp[3])
        raw_seq=tmp[9]
        raw_ins_start=int(tmp[20])
        # left and right seq of refrence and raw reads
        ref_left=ref_seq[raw_ins_start-seq_start-ins_len:raw_ins_start-seq_start]
        ref_right=ref_seq[raw_ins_start-seq_start:raw_ins_start-seq_start+ins_len]
        raw_left=raw_seq[raw_ins_start-ins_len-1:raw_ins_start-1]
        raw_right=raw_seq[raw_ins_start+ins_len-1:raw_ins_start+ins_len*2-1]
        ## detect duplicate position 
        left_ratio=max([similar_equal(ins_seq,ref_left),similar_equal(ins_seq,raw_left)])
        right_ratio=max([similar_equal(ins_seq,ref_right),similar_equal(ins_seq,raw_right)])
        ## if left_ratio==right_ratio: dul_index='right'. consistant with previous program 
        if left_ratio > right_ratio:
            dul_index='left'
        else:
            dul_index='right'
        if dul_index=='left':
            ins_seq_ref=ref_left
            ins_seq_raw=raw_left
            envent=False
            for a in range(len(ins_seq),0,-1):
                ins_ref_left=ins_seq_ref[0:a]
                ins_ref_right=ins_seq_ref[-a:]
                ins_raw_left=ins_seq_raw[0:a]
                ins_raw_right=ins_seq_raw[-a:]
                # attention
                ref_left_seed=ref_seq[raw_ins_start-seq_start-ins_len-a:raw_ins_start-seq_start-ins_len]
                ref_right_seed=ref_seq[raw_ins_start-seq_start:raw_ins_start-seq_start+a]
                seq_left_seed=raw_seq[raw_ins_start-ins_len-1-a:raw_ins_start-ins_len-1]
                seq_right_seed=raw_seq[raw_ins_start+ins_len-1:raw_ins_start+ins_len-1+a]
                if ins_ref_left==ref_right_seed or ins_ref_left==seq_right_seed or ins_raw_left==ref_right_seed or ins_raw_left==seq_right_seed:
                    ins_seed=ins_ref_left
                    seed_len=len(ins_seed)
                    all_ins_microhomo.append(('\t').join([line.strip(),ins_seed,str(seed_len),dul_index,'left']))
                    envent=True
                    break
                elif ins_ref_right==ref_left_seed or ins_ref_right==seq_left_seed or ins_raw_right==ref_left_seed or ins_raw_right==seq_left_seed:
                    ins_seed=ins_ref_right
                    seed_len=len(ins_seed)
                    all_ins_microhomo.append(('\t').join([line.strip(),ins_seed,str(seed_len),dul_index,'right']))
                    envent=True
                    break
        else:
            ins_seq_ref=ref_right
            ins_seq_raw=raw_right
            envent=False
            for a in range(len(ins_seq),0,-1):
                ins_ref_left=ins_seq_ref[0:a]
                ins_ref_right=ins_seq_ref[-a:]
                ins_raw_left=ins_seq_raw[0:a]
                ins_raw_right=ins_seq_raw[-a:]
                # attention
                ref_left_seed=ref_seq[raw_ins_start-seq_start-a:raw_ins_start-seq_start]
                ref_right_seed=ref_seq[raw_ins_start-seq_start+ins_len:raw_ins_start-seq_start+ins_len+a]
                seq_left_seed=raw_seq[raw_ins_start-1-a:raw_ins_start-1]
                seq_right_seed=raw_seq[raw_ins_start+ins_len*2-1:raw_ins_start+ins_len*2-1+a]
                if ins_ref_left==ref_right_seed or ins_ref_left==seq_right_seed or ins_raw_left==ref_right_seed or ins_raw_left==seq_right_seed:
                    ins_seed=ins_ref_left
                    seed_len=len(ins_seed)
                    all_ins_microhomo.append(('\t').join([line.strip(),ins_seed,str(seed_len),dul_index,'left']))
                    envent=True
                    break
                elif ins_ref_right==ref_left_seed or ins_ref_right==seq_left_seed or ins_raw_right==ref_left_seed or ins_raw_right==seq_left_seed:
                    ins_seed=ins_ref_right
                    seed_len=len(ins_seed)
                    all_ins_microhomo.append(('\t').join([line.strip(),ins_seed,str(seed_len),dul_index,'right']))
                    envent=True
                    break            
        if not envent:
            all_ins_microhomo.append(('\t').join([line.strip(),'none','0','none','none']))
f.close()


# %%
prefix=ins_sum.split('/')[-1].split('.')[0]
f=open(path.join(out_put,prefix+'_MH.txt'),'w')
f.write(h1.strip()+'\t'+'ins_MH'+'\t'+'MH_len'+'\t'+'dul_direction'+'\t'+'MH_direction'+'\n')
for i in all_ins_microhomo:
    f.write(i+'\n')
f.close()


# %%



