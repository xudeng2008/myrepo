import os #import os

IDlist=['SRR5660030','SRR5660033','SRR5660044','SRR5660045']

for id in IDlist:
    os.system('fastq-dump -I --split-files '+id)

