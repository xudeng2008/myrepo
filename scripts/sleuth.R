#load library
library(sleuth)

#read in table describing samples and kallisto output
stab <- read.table("sample_table.txt",header=TRUE,stringsAsFactors=FALSE)

#initialize sleuth object
so <- sleuth_prep(stab)

#fit a model comparing the two conditions
so <- sleuth_fit(so, ~condition, 'full')

#fit the reduced model to compare in the likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so, 'reduced', 'full')

#load the dplyr package for data.frame filtering
library(dplyr)

#extract the test results from the sleuth object
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter most significant results (FDR/qval < 0.05)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)

#just show transcript, pval, qval (select by column headers)
sleuth_significant <- dplyr::select(sleuth_significant, target_id,test_stat, pval, qval)

#write data to a file
write.table(sleuth_significant, file="miniProject_Xufang_Deng/DEGs.txt",quote = FALSE,row.names = FALSE)
