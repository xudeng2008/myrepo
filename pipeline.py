
import os #import os
from Bio import Entrez #import Entrez module to retrive data from NCBI
from Bio import SeqIO #import SeqIO module
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#extract sequence data from SRR database
with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list
    for id in IDlist:
        os.system('fastq-dump -I --split-files '+id+" -O miniProject_Xufang_Deng/data") #use fastq-dump to get splited fastq files


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

log = open ("miniProject.log",'w')
log.write('The HCMV genome (EF99921) has '+str(count)+ ' CDS.\n') #write to the log file
log.close()

#Use kallisto to quantify the reads
os.system('mkdir index')
os.system('time kallisto index -i index/index.idx reference_cDNA.fasta') #build an index file

with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list
    for id in IDlist:
        os.system ('time kallisto quant -i index/index.idx -o /miniProject_Xufang_Deng/results/'+id+' -b 30 -t 4 /myrepo/miniProject_Xufang_Deng/'+id+'_1.fastq'+'/miniProject_Xufang_Deng/'+id+'_2.fastq') #run kallisto function to quantify reads

#run sleuth Rscript to calculate differential expression
os.system('Rscript scripts/sleuth.R') #store the significant genes in a DEGs.txt file

with open ("DEGs.txt",'r') as data: # store the output of sleuth to the log file
    with open ("miniProject.log",'a') as log:
        for line in data:
            log.write(line)
log.close()

#Perform mapping using bowtie2
with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
    IDlist = infile.read().splitlines() #make the IDs as list

    os.system('bowtie2-build reference_cDNA.fasta index/HCMV') #build an index file

    for id in IDlist:
        os.system("bowtie2 --quiet -x index/HCMV -1 miniProject_Xufang_Deng/"+id+"_1.fastq"+" -2 miniProject_Xufang_Deng/"+id+"_2.fastq --al-conc al-conc miniProject_Xufang_Deng/results/"+id+"mapped.fq") #map reads to a reference


#define a function to count the number of reads
def count_read(infile):
    count = 0
    for line in infile:
        count += 1
    read_count = count/4 #each read has 4 lines
    return read_count


with open ("miniProject.log",'a') as log:
    with open ("SRR_accession_IDs.txt",'r') as infile: #read in the accession IDs
        IDlist = infile.read().splitlines() #make the IDs as list
        donorIDs = ['1','1','3','3']
        timePoints = ['2','6','2','6']
        i = 0
        for id in IDlist:
            infile1 = open('miniProject_Xufang_Deng/'+id+'_1.fastq','r')
            total_raw = count_read(infile1)
            
            infile2 = open('miniProject_Xufang_Deng/results/'+id+'mapped.1.fq','r')
            total_mapped = count_read(infile2)
            log.write("Donor "+donorIDs[i]+" at "+timePoints[i]+" dpi had "+str(total_raw)+" read pairs before Bowtie2 filtering and "+str(total_mapped)+" read after.\n\n")
            i+=1
log.close()

#Use SPAdes to assemble mapped reads
os.system('cat miniProject_Xufang_Deng/results/*mapped.1.fq > merged_1.fq') #merge all 4 forward file to  a single file
os.system('cat miniProject_Xufang_Deng/results/*mapped.2.fq > merged_2.fq') #merge all 4 reverse file to  a single file
os.system('spades -k 55,77,99,127 -t 2 --only-assembler -1 miniProject_Xufang_Deng/results/merged_1.fq -2 miniProject_Xufang_Deng/results/merged_2.fq  -o HCMV_assembly') #use SPAdes to assembly 

with open ("miniProject.log",'a') as log:
    log.write('spades -k 55,77,99,127 -t 2 --only-assembler -1 miniProject_Xufang_Deng/results/merged_1.fq -2 miniProject_Xufang_Deng/results/merged_2.fq -o HCMV_assembly\n\n')
log.close()

#count the number of contig with a length greater than 1000 nt
#extract these contigs
with open("HCMV_assembly/contigs.fasta",'r') as infile:
    with open ("HCMV_assembly/contigs_gt_1000.fasta",'w') as outfile:
        totalLen = 0
        count = 0
        for record in SeqIO.parse(infile,'fasta'):
            ctgLen = len(record.seq)
            if ctgLen > 1000:
                count += 1
                totalLen += ctgLen
                SeqIO.write(record,outfile, 'fasta') #put the contigs longer than 1000bp in a file
with open ("miniProject.log",'a') as log:
    log.write("There are "+str(count)+" contigs > 1000 bp in the assembly.\n\n")
    log.write("There are "+str(totalLen)+" bp in the assembly.\n\n")
log.close()
#concatenate all of the contigs > 1000 bp with a linker
#build a linker string
linker = ''
while len(linker)<50:
    linker += 'N'

#concatenate all of the contigs > 1000 bp
with open('HCMV_assembly/contigs_gt_1000.fasta','r') as infile:
    with open ('HCMV_assembly/ligated_contigs.fasta','w') as outfile:
        superString = ''
        for record in SeqIO.parse(infile, 'fasta'):
            currCtg = str(record.seq)+linker #add the linker
            superString += currCtg
    outfile.write(">superString\n"+superString[:len(superString)-50]) #remove the last linker from the string

outfile.close()

#Blast the assembled sequence in NCBI
query_seq = open('ligated_contigs.fasta').read()
result_handle = NCBIWWW.qblast("blastn", "nt", query_seq, entrez_query = "Herpesviridae", hitlist_size= 10)

with open ("miniProject.log",'a') as log:
    log.write("seq_title\talign_len\tnumber_HSPs\ttopHSP_ident\ttopHSP_gaps\ttopHSP_bits\ttopHSP_expect")

blast_records = NCBIXML.parse(result_handle)
blast_records = list(blast_records)
for record in blast_records:
    for alignment in record.alignments:
        seq_title= alignment.title
        align_len= str(alignment.length)
        number_HSPs= str(len(alignment.hsps))
        topHSP_ident= str(alignment.hsps[0].identities)
        topHSP_gaps= str(alignment.hsps[0].gaps)
        topHSP_bits= str(alignment.hsps[0].bits)
        topHSP_expect= str(alignment.hsps[0].expect)
        log.write(seq_title+"\ta"+lign_len+"\t"+number_HSPs+"\t"+topHSP_ident+"\t"+topHSP_gaps+"\t"+topHSP_bits+"\t"+topHSP_expect)

log.close()
