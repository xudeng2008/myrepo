from Bio import SeqIO

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

