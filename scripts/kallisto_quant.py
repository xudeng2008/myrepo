import os
#Use kallisto to quantify the reads
os.system('time kallisto index -i index/index.idx reference_cDNA.fasta') #build an index file

with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list
    for id in IDlist:
        os.system ('time kallisto quant -i index/index.idx -o /miniProject_Xufang_Deng/results/'+id+' -b 30 -t 4 /myrepo/miniProject_Xufang_Deng/'+id+'_1.fastq'+'/miniProject_Xufang_Deng/'+id+'_2.fastq') #run kallisto function to quantify reads

