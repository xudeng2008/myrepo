#import Entrez module to retrive data from NCBI
from Bio import Entrez
from Bio import SeqIO
#Fill in the Entrez.email field
Entrez.email = "xudeng@luc.edu.com"


#retrive data by searching the Nucleotide database with term
handle = Entrez.efetch(db="nucleotide", id='EF999921', rettype="gb", retmode="text")

with open("reference_cDNA.fasta",'w') as outfile:
    for record in SeqIO.parse(handle,"genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == 'CDS':
                    pdtlist = feature.qualifiers["product"]
                    outfile.write(">"+pdtlist[0]+'\n')
                    CDSlist = feature.location.extract(record).seq
                    outfile.write(str(CDSlist)+'\n')
    outfile.close()
