import os
IDlist = ['SRR5660030','SRR5660033','SRR5660044','SRR5660045']
for id in IDlist:
        os.system ('time kallisto quant -i index.idx -o /homes/xudeng/myrepo/miniProject_Xufang_Deng/results/'+id+' -b 30 -t 4 /homes/xudeng/myrepo/miniProject_Xufang_Deng/'+id+'_1.fastq'+' /homes/xudeng/myrepo/miniProject_Xufang_Deng/'+id+'_2.fastq')

