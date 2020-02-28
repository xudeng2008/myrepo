from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#Blast the assembled sequence in NCBI
query_seq = open('HCMV_assembly_test/ligated_contigs.fasta').read() # read in the assembled contigs
result_handle = NCBIWWW.qblast("blastn", "nr", query_seq, entrez_query = "Herpesviridae", hitlist_size= 10) # run blastn in nt database

with open ("miniProject_test.log",'a') as log: #read in the log file
    log.write("seq_title\talign_len\tnumber_HSPs\ttopHSP_ident\ttopHSP_gaps\ttopHSP_bits\ttopHSP_expect\n") # head of the log file of blast report

    blast_records = NCBIXML.parse(result_handle) #Use NCBIXML to parse the blast output
    blast_records = list(blast_records) #convert the output as a list
    for record in blast_records: #loop the list of blast results
        for alignment in record.alignments: #loop each alignment
            seq_title= alignment.title #get the subject sequence title
            align_len= str(alignment.length) # alignment length
            number_HSPs= str(len(alignment.hsps)) # number of high-scoring segments in the alignment
            topHSP_ident= str(alignment.hsps[0].identities) # the identity between query sequence and subjected segment
            topHSP_gaps= str(alignment.hsps[0].gaps) # find the gaps of the alignment
            topHSP_bits= str(alignment.hsps[0].bits) # bit-score of the alignment
            topHSP_expect= str(alignment.hsps[0].expect) # expect value of the alignment
            # write these output in log file
            log.write(seq_title+"\t"+align_len+"\t"+number_HSPs+"\t"+topHSP_ident+"\t"+topHSP_gaps+"\t"+topHSP_bits+"\t"+topHSP_expect+'\n')
#close the log file
log.close()
