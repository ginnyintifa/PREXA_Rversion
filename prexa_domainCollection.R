## 2 gethering the domain indices on the same set of proteins

library(dplyr)
library(magrittr)
library(data.table)



domain_prots=fread("/Users/ginny/tcga_june_gastric_cancer/domain/protid.tsv",header=F)
domain_name=fread("/Users/ginny/tcga_june_gastric_cancer/domain/domain_name.tsv",header=F)
domain_start=fread("/Users/ginny/tcga_june_gastric_cancer/domain/env_start.tsv",header=F)
domain_end=fread("/Users/ginny/tcga_june_gastric_cancer/domain/env_end.tsv",header=F)
domain_type=fread("/Users/ginny/tcga_june_gastric_cancer/domain/domain_type.tsv",header=F)

###what kind of output do I want?

setwd("/Users/ginny/tcga_june_gastric_cancer")

domain_df=data.frame(domain_prots, domain_start, domain_end, domain_name, domain_type)

colnames(domain_df)=c("ID","start_position","end_position","domain_name","domain_type")
write.table(domain_df, "domain_df.tsv", quote=F, row.names = F, sep="\t")




protein_id_dic=unique(data.frame(sorted_add_onlyID$protein_ID, sorted_add_onlyID$ID))
colnames(protein_id_dic)=c("protein_ID","ID")

###find the full id for these seqs\
which_is_in=which(domain_df$protein_ID %in% unique(sorted_add_onlyID$ID))

domain_df_sp=domain_df[which_is_in, ]

domain_df_sp$protein_ID=protein_id_dic$protein_ID[match(domain_df_sp$ID,protein_id_dic$ID)]

domain_frame=domain_df_sp[,c(6,1,2,3,4,5)]

domain_sorted=domain_frame[order(domain_frame$protein_ID),]

write.table(domain_sorted, "domain_sorted_df.tsv", quote=F, row.names = F, sep="\t")
