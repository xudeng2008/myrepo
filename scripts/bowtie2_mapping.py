import os

#Perform mapping using bowtie2
with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list

    os.system('bowtie2-build '+path+'reference_genome.fasta index_test/HCMV_test') #build an index file

    for id in IDlist:
        os.system("bowtie2 --quiet -x index_test/HCMV_test -1 "+path+id+"_1.fastq"+" -2 "+path+id+"_2.fastq --al-conc "+path+id+"mapped.fq") #map reads to a reference


