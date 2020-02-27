from Bio import Entrez #import Entrez module to retrive data from NCBI
from Bio import SeqIO #import SeqIO module

#Fill in the Entrez.email field
Entrez.email = "xudeng@luc.edu"
#retrive data by searching the Nucleotide database with term
handle = Entrez.efetch(db="nucleotide", id='EF999921', rettype="gb", retmode="text")
with open("reference_cDNA.fasta",'w') as outfile: #create an out file
    count = 0 #store the number of CDS
    for record in SeqIO.parse(handle,"genbank"): #parse the file
        if record.features: #if record has features
            for feature in record.features: #loop the features
                if feature.type == 'CDS': #find the CDS feature
                    pdtlist = feature.qualifiers["product"] #get the product name
                    outfile.write(">"+pdtlist[0]+'\n') #store the name
                    CDSlist = feature.location.extract(record).seq #get the CDS sequence
                    outfile.write(str(CDSlist)+'\n')
                    count+=1
    outfile.close()
