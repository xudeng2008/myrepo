from Bio import SeqIO

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
