### getting ptm first data

## gethering PTM info from PTMtopographer


library(data.table)
library(dplyr)
library(magrittr)


setwd("/Users/ginny/tcga_june_gastric_cancer")


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




