from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

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
        log.write(seq_title\talign_len\tnumber_HSPs\ttopHSP_ident\ttopHSP_gaps\ttopHSP_bits\ttopHSP_expect\n)

log.close()
