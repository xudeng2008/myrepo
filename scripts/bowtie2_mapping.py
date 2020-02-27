import os

#Perform mapping using bowtie2
with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list

    os.system('bowtie2-build reference_cDNA.fasta index/HCMV') #build an index file

    for id in IDlist:
        os.system("bowtie2 --quiet -x index/HCMV -1 miniProject_Xufang_Deng/"+id+"_1.fastq"+" -2 miniProject_Xufang_Deng/"+id+"_2.fastq --al-conc al-conc miniProject_Xufang_Deng/results/"+id+"mapped.fq") #map reads to a reference

