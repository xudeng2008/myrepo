import os #import os

#extract sequence data from SRR database
with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list
    for id in IDlist:
        os.system('fastq-dump -I --split-files '+id+" -O miniProject_Xufang_Deng/data") #use fastq-dump to get splited fastq files

