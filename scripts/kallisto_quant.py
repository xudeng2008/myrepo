import os

os.system('mkdir index_test')
os.system('time kallisto index -i index_test/index.idx '+path+'reference_cDNA.fasta') #build an index file

with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list
    for id in IDlist:
        os.system ('time kallisto quant -i index_test/index.idx -o '+path+id+' -b 30 -t 4 '+path+id+'_1.fastq '+path+id+'_2.fastq') #run kallisto function to quantify reads
