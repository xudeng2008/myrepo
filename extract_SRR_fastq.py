import os #import os
with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
  IDlist = infile.read().splitlines() #make the IDs as list
for id in IDlist:
    os.system('fastq-dump -I --split-files '+id) #use fastq-dump to get splited fastq files

