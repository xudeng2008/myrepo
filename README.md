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
5. fastq-dump
6. And all their dependencies

## Scripts:
**1. extract_reference_cDNA.py** Extract the reference cDNAs from NCBI with given accession IDs. This will generate a fasta file to store the output. The output file will be used by kallisto to build an index.
**2. extract_reference_genome.py** Extract the reference genome from NCBI with given accession IDs. This will generate a fasta file to store the output. The output file will be used by Bowtie2 to build an index.
**3. extract_SRR_fastq.py** Extract sequencing data from the SRA database with given accession numbers.
**4. kallisto_quant.py** Use kallisto to quantify the read counts. This script includes building an index.
**5. sleuth.R** An R package uses kallisto output to perform statistical tests for differential expression. This will generate a DEGs.txt file to store the differential expressed genes (q<0.05).
**6. bowtie2_mapping.py** Use bowtie2 to map reads to reference sequences. This script includes building an index.
**7. SPAdes_assembly.py** Use SPAdes program to assemble mapped reads.
**8. extract_contigs_gt1000.py** This script uses a contigs.fasta file of SPAdes output to extract assembled contigs with 1000 bp or longer. A fasta file is generated. 
**9. ligate_contigs.py** This scripts uses the contigs output file (#7) to assemble a large contig.
**10. blast_assembly.py** Use the assembled large contig (#8) to blast in Genbank.

## Python modules:
1. os
2. Biopython: SeqIO, Entrez, NCBIWWW, NCBIXML

## Supplemental files:
**1. SRR_accession_IDs.txt** A file contains accession numbers
**2. sample_table.txt** A file contains a table for sleuth package doing statitic analysis
