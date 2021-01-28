# -*- encoding: utf-8 -*-
'''
@File    :   total_reads.py
@Time    :   2020/07/21 11:03:23
@Author  :   Binbin Wang 
@Version :   1.0
@Contact :   wangbinbintj@gmail.com
'''
#%%
# here put the import lib
import os
from os import path
import sys
#%%
filt_clones_files=sys.argv[1]
out_put=sys.argv[2]

if not path.exists(out_put):
    os.makedirs(out_put)

# %%
# append total reads info for new samples
total_reads_sum={}
for i in os.listdir(filt_clones_files):
    if not i.startswith('.') and i.endswith('filt_clones.txt'):
        print(i)
        f=path.join(filt_clones_files,i)
        tmp_name=i[0:6]
        if tmp_name not in total_reads_sum:
            tmp_reads=len(open(f).readlines())-1
            total_reads_sum[tmp_name]=str(tmp_reads)
            
#%%
f=open(out_put+'/total_reads.txt','w')
for i in total_reads_sum:
    f.write(('\t').join([i,total_reads_sum[i]])+'\n')
f.close()
#%%
