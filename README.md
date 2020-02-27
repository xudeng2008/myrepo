# Goal: This pipeline scripts is designed for analysing sequencing data, assembling a genome, and blasting a genome in Genbank database. 

## Environment Requirements:
1. Linux
2. Python 3.6
3. R

## Programs:
1. Bowtie2
2. Kallisto
3. Sleuth
4. SPAdes
5. And all their dependencies

## Scripts:
1. **extract_reference_cDNA.py** Extract the reference cDNAs from NCBI with given accession IDs. This will generate a fasta file to store the output.
2. **extract_reference_genome.py** Extract the reference genome from NCBI.
3. **kallisto_quant.py** Use kallisto to quantify the read counts. This script includes building an index.
4. **bowtie2_mapping.py** Use bowtie2 to map reads to reference sequences. This script includes building an index.
5. **SPAdes_assembly.py** Use SPAdes program to assemble mapped reads.
6. **extract_contigs_gt1000.py** Extract assembled contigs with 1000 bp or longer.
7. **sleuth.R** An R package uses kallisto output to perform statistical tests for differential expression
8. **ligate_contigs.py** Ligate contigs to assemble a large contig.
9. **blast_assembly.py** Use the assembled large contig to blast in Genbank.
10. **extract_SRR_fastq.py** Extract sequencing data from the SRA database with given accession numbers.
