from Bio import Entrez #import Entrez module to retrive data from NCBI
from Bio import SeqIO #import SeqIO module

#Fill in the Entrez.email field
Entrez.email = "youremail@luc.edu"

#retrive data by searching the Nucleotide database with term
handle = Entrez.efetch(db="nucleotide", id='EF999921', rettype="gb", retmode="text")
with open("miniProject_Xufang_Deng/reference_cDNA.fasta",'w') as outfile: #create an out file
#    record = SeqIO.parse(handle,"genbank")
#    record.features = [f for f in record.features if f.type == "CDS"]
#    SeqIO.write(record, outfile, "genbank")
#    outfile.close()
    for record in SeqIO.parse(handle,"genbank"): #parse the file
        if record.features: #if record has features
            for feature in record.features: #loop the features
                if feature.type == 'CDS': #find the CDS feature
                    outfile.write(">"+str(feature.qualifiers["protein_id"]).replace("[","").replace("]","").replace("'","")+'\n') #store the name of the CDS
                    CDS = feature.location.extract(record).seq #get the CDS sequence
                    outfile.write(str(CDS)+'\n')
    outfile.close()
