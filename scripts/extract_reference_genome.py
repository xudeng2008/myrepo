from Bio import SeqIO
from Bio import Entrez
email = "youremail@luc.edu"
with open ("miniProject_Xufang_Deng/reference_genome.fasta",'w') as outfile: 
    #Fill in the Entrez.email field
    Entrez.email = email
    #retrive data by searching the Nucleotide database with term
    handle = Entrez.efetch(db="nucleotide", id='EF999921', rettype="fasta")
    record = SeqIO.read(handle, "fasta")
    outfile.write(">"+str(record.description)+'\n'+str(record.seq)+'\n')
    outfile.close()
