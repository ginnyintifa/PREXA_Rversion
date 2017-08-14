
### Find hyper mutated and patients do not have as many, exclude them 


### the following, I'm looking at the sum of mutations each patient has and
### exclude hypermutated people.the criterion is derived from tcga paper
### and derived by hw



#### major change 

### I think in this section I'm setting upper and lower limit of the patients selection, exclude hypermutated patients and exclude patients with too few mutations.


#setwd("/Users/ginny/tcga_june_gastric_cancer/TCGA_STAD")


### the premise is stad whole is on the same order of stad clinical 



get_selected_cox_matrix = function (stad_whole, stad_clinical, h_cutoff, l_cutoff, cancer_type)
{
  
  stad_matrix = stad_whole[,-1]  #### the first column is patient barcode
  
  sum_mut_per_patient = apply(stad_matrix,1,sum)
  
  
  #cutoff = 3500 ### derived by hw
  
  sel_which = which(sum_mut_per_patient <= h_cutoff & sum_mut_per_patient >= l_cutoff)
  
  sel_patient = stad_clinical[sel_which,]
  
  ### get the order of patients right in next section
  
  write.table(sel_patient, paste(cancer_type,"sel_clinical.tsv", sep = "_"),
              sep = "\t", quote = F, row.names = F)
  
  

  ### arrange by sur_ind and sur_time
  stad_sel_sort_clinical <- sel_patient %>%
    arrange(survival_ind,desc(survival_time))
   
  write.table(stad_sel_sort_clinical, paste(cancer_type,"sel_sort_clinical.tsv",sep = "_"),
              sep = "\t",quote = F, row.names = F)
  
  
  
  cat (table(sel_patient$survival_ind))
  
  
  sel_matrix = stad_matrix[sel_which,]
  sel_matrix_patient = data.frame(sel_patient$bcr_patient_barcode,sel_matrix)
  colnames(sel_matrix_patient) = colnames(stad_whole)
  write.table(sel_matrix_patient, paste(cancer_type,"sel_cox_table.tsv", sep = "_"),sep = "\t", quote = F, row.names = F)
  
  
  t_sel_matrix = t(sel_matrix)
  t_sel_df = data.frame(colnames(stad_whole)[-1],t_sel_matrix)  
  
  colnames(t_sel_df) = c("locus",sel_patient$bcr_patient_barcode)
  
  write.table(t_sel_df, paste(cancer_type, "sel_cox_table_transpose.tsv", sep="_"),sep = "\t", quote = F, row.names = F)
  
  
}



#stad_whole = fread("stad_full_cox_table.tsv")
#stad_clinical = fread("stad_clinical.tsv",stringsAsFactors = F)



##table(apply(stad_whole[,-1],1,sum))

### visualize the total number of mutations by histogram
##hist(apply(stad_whole[,-1],1,sum),main = "sum of mutations each person has", breaks = 50)


### set cutoff of hypermutaion patients at 3500. derived from the tcga gastric paper
#get_selected_cox_matrix(stad_whole,  stad_clinical, 3500, 1, "stad")

### will also produce a table of selected patient info 



