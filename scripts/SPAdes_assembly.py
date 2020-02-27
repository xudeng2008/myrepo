import os

os.system('cat miniProject_Xufang_Deng/results/*mapped.1.fq > merged_1.fq') #merge all 4 forward file to  a single file
os.system('cat miniProject_Xufang_Deng/results/*mapped.2.fq > merged_2.fq') #merge all 4 reverse file to  a single file
os.system('spades -k 55,77,99,127 -t 2 --only-assembler -1 miniProject_Xufang_Deng/results/merged_1.fq -2 miniProject_Xufang_Deng/results/merged_2.fq  -o HCMV_assembly') #use SPAdes to assembly 

