import os

#Use SPAdes to assemble mapped reads
os.system('cat '+path+'*mapped.1.fq > '+path+'merged_1.fq') #merge all 4 forward file to  a single file
os.system('cat '+path+'*mapped.2.fq > '+path+'merged_2.fq') #merge all 4 reverse file to  a single file
os.system('spades -k 55,77,99,127 -t 4 --only-assembler -1 '+path+'merged_1.fq -2 '+path+'merged_2.fq  -o HCMV_assembly_test') #use SPAdes to assembly 


