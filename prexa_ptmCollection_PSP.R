### getting ptm first data

## gathering PTM info from PSP positive sites only

### Sunday I will focus on this

### I will produce the BRCA data on Sunday


library(data.table)
library(dplyr)
library(magrittr)


#setwd("/Users/ginny/PREXA_Rproject/")


## I plan to re map the PSP data to the sp proteins 
## use the korean mappign program to do this 


### preprocessing of the pep data from PSP, only the middle letter is considered as 
### modified

 
## 2 criteria

##ok my concern is 

### if I map to all I need to consider the tail and head problem
### i can use the map to one strategy to map to only the protein mentioned


### so for now the idea is to map only the 15 mers to proteins at the loss of missing 
### head and tail

### the first thing to do is to make sure the peptide file is correct, only the middle
### letter of the 15mers are lower letters


### ok the thing is to get the exact mapping from the PTM topographer software
### and this version is for exact mapping

### get the exact mapping version of PTMtopographer

### new PSP data?



## I think the best way to do this is to look at only sites from PhosphositePlus and
## the fasta seq file where PTM sites are indicated in lowercase letters.


library(data.table)
library(magrittr)
library(dplyr)



## gethering phosphorylation first
phospho_fasta = readLines("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/Phosphosite_PTM_seq.fasta")

## get the humamn proteins


phos_id_ind = which(substr(phospho_fasta,1,1)==">")
phos_id = phospho_fasta[phos_id_ind]
## get the corresponding sequence
seq_id = rep(list(c()),length(phos_id))
for (i in 1:(length(seq_id)-1))
{
  seq_id[[i]]= c((phos_id_ind[i]+1):(phos_id_ind[i+1]-1))
  
}
seq_id[[length(seq_id)]] = c((phos_id_ind[length(seq_id)]+1):length(phospho_fasta))
## concatenate
seq = unlist(sapply(1:length(phos_id), function(x) 
  paste0(phospho_fasta[seq_id[[x]]],collapse = "")))

phos_species = sapply(1:length(phos_id), function(x) 
  unlist(strsplit(phos_id[x],split = "|",fixed = T))[3])
phos_id_human_ind = which(phos_species == "human")

phos_id_human = phos_id[phos_id_human_ind]
phos_seq_human = seq[phos_id_human_ind]

sum_s = 0

for (i in 1:length(phos_seq_human))
{
  sum_s = sum_s + sum(sapply(1:nchar(phos_seq_human[i]), function(x)
    substr(phos_seq_human[i],x,x) == "s"))
  

}

sum_t = 0

for (i in 1:length(phos_seq_human))
{
  sum_t = sum_t + sum(sapply(1:nchar(phos_seq_human[i]), function(x)
    substr(phos_seq_human[i],x,x) == "t"))
  
  
}

sum_y = 0

for (i in 1:length(phos_seq_human))
{
  sum_y = sum_y + sum(sapply(1:nchar(phos_seq_human[i]), function(x)
    substr(phos_seq_human[i],x,x) == "y"))
  
  
}



###########

## write out the file of phosphorylation fasta

phos_file = rep("",2*length(phos_seq_human))
for (i in 1:length(phos_seq_human))
{
  phos_file[2*i-1] = phos_id_human[i]
  phos_file[2*i] = phos_seq_human[i]
  
}


write.table(phos_file, "phosphorylation_site_fasta.tsv", 
            quote = F, row.names = F)


### on top of this, I will identify the position of each PTM site and generate
### the classic PTM file for korean gastric mapping

library(data.table)
library(magrittr)
library(dplyr)


### extract position and prot ID

### in the format of PTM_aggregate_PSP

phos_data = readLines("/Users/ginny/Korean_gastric/phosphorylation_site_fasta.tsv")

phos_id = phos_data[c(T,F)]
phos_seq = phos_data[c(F,T)]


phos_id_uni = sapply(1:length(phos_id), 
                     function(x) strsplit(phos_id[x], split = "|",fixed = T)[[1]][4])

all_info = c("","","")

for (i in 1:length(phos_seq))
{
  s_seq = sapply(1:nchar(phos_seq[i]), function(x)
    substr(phos_seq[i],x,x) == "s")
  
  t_seq = sapply(1:nchar(phos_seq[i]), function(x)
    substr(phos_seq[i],x,x) == "t")
  
  y_seq = sapply(1:nchar(phos_seq[i]), function(x)
    substr(phos_seq[i],x,x) == "y")
  
  all_pos = c(1:nchar(phos_seq[i]))
  
  s_pos = all_pos[s_seq]
  t_pos = all_pos[t_seq]
  y_pos = all_pos[y_seq]
  
  s_info = 
    cbind(rep(phos_id_uni[i],length(s_pos)), s_pos, rep("ps", length(s_pos)))
  t_info = 
    cbind(rep(phos_id_uni[i],length(t_pos)), t_pos, rep("pt", length(t_pos)))
  y_info = 
    cbind(rep(phos_id_uni[i],length(y_pos)), y_pos, rep("py", length(y_pos)))
  
  
  all_info = rbind(all_info, s_info, t_info, y_info)
  
  
  
  
}

phospho_info = all_info[-1,]
phospho_df = data.frame(ptm_ID = phospho_info[,1], 
                        ptm_position = phospho_info[,2],
                        ptm_type = phospho_info[,3])



write.table(phospho_df, "ptm_aggregate_PSP.tsv", quote = F, row.names = F, sep = "\t")


##################################################################
# I can use the following to process any type of ptm data 


####

## I will work on collect other PTMs today because others are not in fasta
## sequences


### extract ubiquitination



parse_PSP = function(ubi_psp, mod_site, type)
{
  #ubi_psp = sumoy_psp
  #mod_site = "K"
  #type = "sumoy"
  
  ubi_site <- ubi_psp %>% 
    filter(ORGANISM == "human") %>% 
    select(ACC_ID, MOD_RSD)
  
  ### make sure only select K not any other
  
  ubi_sel = ubi_site[which(substr(ubi_site$MOD_RSD,1,1)== mod_site),]
  
  it = which(substr(ubi_site$MOD_RSD,1,1)!= mod_site)
  
  ubi_residue = sapply(1:nrow(ubi_sel), 
                       function(x) strsplit(ubi_sel$MOD_RSD[x],"-")[[1]][1])
  
  
  ### use regular experssion to get rid of the letters
  
  
  
  ptm_position = as.numeric(gsub("[A-Z]","",ubi_residue))
  
  ptm_type = rep(type, length(ptm_position))
  
  
  ubi_df = data.frame(ptm_ID = ubi_sel$ACC_ID, ptm_position, ptm_type)
  
  
  write.table(ubi_df, paste0(type,"_PSP.tsv"), quote = F, row.names = F, sep = "\t")
  
  return (ubi_df)
}

### filter out other species and select only acc id and residue position




ubi_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/Ubiquitination_site_dataset")
ubi_df = parse_PSP(ubi_psp, "K", "ubi")


phospho_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/Phosphorylation_site_dataset")
ps_df = parse_PSP(phospho_psp, "S", "ps")
pt_df = parse_PSP(phospho_psp, "T", "pt")
py_df = parse_PSP(phospho_psp, "Y", "py")

acety_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/Acetylation_site_dataset")
acety_df = parse_PSP(acety_psp, "K", "acety")

methy_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/Methylation_site_dataset")
methy_k_df = parse_PSP(methy_psp, "K", "methy_k")
methy_r_df = parse_PSP(methy_psp, "R", "methy_r")


sumoy_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/Sumoylation_site_dataset")
sumoy_df = parse_PSP(sumoy_psp, "K", "sumoy")

glcnac_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/O-GlcNAc_site_dataset")
glcnac_s_df = parse_PSP(glcnac_psp, "S", "glcnac_s")
glcnac_t_df = parse_PSP(glcnac_psp, "T", "glcnac_t")


galnac_psp = fread("/Users/ginny/PTMtopographer_2017_08/PSP_Downloads_FEB_16_2017/O-GalNAc_site_dataset")
galnac_s_df = parse_PSP(galnac_psp, "S", "galnac_s")
galnac_t_df = parse_PSP(galnac_psp, "T", "galcnac_t")


#### combine

all_ptm_df = rbind(ubi_df, ps_df, pt_df, py_df, 
                   acety_df, methy_r_df,methy_k_df,
                   sumoy_df, glcnac_t_df, glcnac_s_df,
                   galnac_t_df, galnac_s_df)


write.table(all_ptm_df, "PSP_ptm_aggregate.tsv", quote =F, row.names = F, sep = "\t")





setwd("/Users/ginny/PREXA_Rproject")

phospho_pep = fread("phospho_pep.tsv", header = F)

center_pep <- phospho_pep%>% filter(nchar(V1)==15)

C_center_pep = toupper(center_pep$V1)

substr(C_center_pep, 8,8) = substr(center_pep$V1,8,8)

### maybe just simply map all to sp, discuss with hw about this
### I think the better way is to only look at those matched protein IDs




###get the sites with in 5% der from both and combine them together

get_df=function(ps_rf, ps_svm,ptm_type)
{
  ps_rf_der=apply(cbind(ps_rf$global_dfdr,ps_rf$pspc_dfdr),1,min)
  
  select_five_percent_rf =which(ps_rf_der<=0.05)
  
  ps_rf_select=cbind(ps_rf$protein_ID[select_five_percent_rf],
                     ps_rf$position[select_five_percent_rf],
                     ps_rf$psp_state[select_five_percent_rf])
  
  ps_svm_der=apply(cbind(ps_svm$global_dfdr,ps_svm$pspc_dfdr),1,min)
  
  select_five_percent_svm =which(ps_svm_der<=0.05)
  ps_svm_select=cbind(ps_svm$protein_ID[select_five_percent_svm],
                      ps_svm$position[select_five_percent_svm],
                      ps_svm$psp_state[select_five_percent_svm])
  
  ps_select=unique(rbind(ps_rf_select, ps_svm_select))
  
  ps_select_df=as.data.frame(cbind(ps_select,rep(ptm_type,nrow(ps_select) )))
  
  colnames(ps_select_df)=c("protein_ID","position","psp_state","ptm_type")
  
  return(ps_select_df)
  
}


ps_rf=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/ps/ps_rf_sp_combine.tsv")
ps_svm=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/ps/ps_svm_sp_combine.tsv")

ps_select_df=get_df(ps_rf, ps_svm,"ps")



pt_rf=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/pt/pt_rf_sp_combine.tsv")
pt_svm=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/pt/pt_svm_sp_combine.tsv")

pt_select_df=get_df(pt_rf, pt_svm,"pt")


py_rf=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/py/py_rf_sp_combine.tsv")
py_svm=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/py/py_svm_sp_combine.tsv")
py_select_df=get_df(py_rf, py_svm,"py")

ubi_rf=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/ubi/ubi_rf_sp_combine.tsv")
ubi_svm=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/cv_phospho_ubi/ubi/ubi_svm_sp_combine.tsv")

ubi_select_df=get_df(ubi_rf, ubi_svm,"ubik")

#########get the selection directly for the other modifications



get_new_df=function(acety_select,ptm_type)
{
  acety_select_df=as.data.frame(cbind(acety_select$protein_ID,
                                      acety_select$position,
                                      acety_select$psp_state,
                                      rep(ptm_type,nrow(acety_select))))
  
  colnames(acety_select_df)=c("protein_ID","position","psp_state","ptm_type")
  
  return(acety_select_df)
  
  
}
acety_select=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/acety_pred/acety_sp_prediction.tsv")

acety_select_df=get_new_df(acety_select, "acety_k")


methy_k_select=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/methy_k_pred/methy_k_sp_prediction.tsv")

methy_k_select_df=get_new_df(methy_k_select, "methy_k")


methy_r_select=fread("/Users/ginny/PTMtopographer_DEC/bin/bin_folder_not_in_git/methy_r_pred/methy_r_sp_prediction.tsv")

methy_r_select_df=get_new_df(methy_r_select, "methy_r")

###then I will merge these ptm together to get my ptm file again
ptm_data=as.data.frame(rbind(ps_select_df,pt_select_df,py_select_df,ubi_select_df,
                             acety_select_df,methy_k_select_df,methy_r_select_df))

sorted_ptm_data=ptm_data[order(ptm_data$protein_ID,ptm_data$position),]


get_id_pos=cbind(sorted_ptm_data$protein_ID,sorted_ptm_data$position)

write.table(sorted_ptm_data, "ptm_data_original.tsv", quote=F, row.names = F, sep="\t")




###now I want to create another column recording the ptm type for each site in each position

keeps=c("protein_ID","position","ptm_type")

get_merge=sorted_ptm_data[keeps]


add_together=function(alist)
{
  return(paste(alist, collapse = "+"))
}


tryit=aggregate(ptm_type~protein_ID+position, data=get_merge,add_together )

sorted_aggregate=tryit[order(tryit$protein_ID),]


onlyID=sapply(sorted_aggregate$protein_ID, function(x) strsplit(as.character(x), split="|",fixed=T)[[1]][2])

sorted_add_onlyID=data.frame(sorted_aggregate$protein_ID,onlyID, sorted_aggregate$position,sorted_aggregate$ptm_type)

colnames(sorted_add_onlyID)=c("protein_ID","ID","position","ptm_type")

write.table(sorted_add_onlyID, "ptm_data_aggregate.tsv", quote=F, row.names = F, sep="\t")




