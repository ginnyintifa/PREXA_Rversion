### get the whole mutation loci matrix for the later use of cox regression
### to get rid of patients who are hyper mutated

get_this_column = function(i, all_locus, this_locus_patient, sel_clinical)
{
  
  
  this_column = rep(0, length(sel_clinical$bcr_patient_barcode))
  get_patient <- this_locus_patient %>% filter(locus_name == all_locus[i])%>% select(sampleid)
  find_patient = which(sel_clinical$bcr_patient_barcode %in% get_patient$sampleid)
  this_column[find_patient]=1
  
  return (this_column)
  
}




get_locus_patient_matrix = function (sel_loci, cancer_type, sel_clinical)
{
  
 # cancer_type = "stad"
  
  #this_locus_patient <- sel_loci %>% mutate(locus_name = paste(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos, mutationaa,aapos,protlen,sep="-")) %>% ungroup() %>% select(locus_name, sampleid)
  
  
  locus_name = paste(sel_loci$genename, sel_loci$uni_swiss, sel_loci$uni_trembl,sel_loci$chromosome, sel_loci$mutationtype, 
                     sel_loci$mutationpos, sel_loci$mutationaa,sel_loci$aapos,
                     sel_loci$protlen,sep="-")
  
  this_locus_patient = data.frame(locus_name, sampleid = sel_loci$sampleid)
  #this_matrix = matrix(0, length(right_patient_info$bcr_patient_barcode), length(all_locus))
  
  all_locus = unique(locus_name)
  the_matrix = sapply(1:length(all_locus), function(i){
    it = get_this_column(i, all_locus, this_locus_patient,sel_clinical)
    if(i%%1000 ==0)
    {
      cat(i)
      cat("\n")
    }
  
    return(it)
    
  } )
  
  
  
  ### figure out why the following step is necessary ### because patient is a subset?
  
  the_matrix_keep_mutation = which(apply(the_matrix,2,sum)>=10)
  
  keep_matrix = the_matrix[,the_matrix_keep_mutation]
  keep_locus = all_locus[the_matrix_keep_mutation]
  
  
  keep_df = data.frame(sel_clinical$bcr_patient_barcode,keep_matrix)
  colnames(keep_df) = c("patient_id",keep_locus)
  
  #cancer_type = "stad"
  write.table(keep_df, paste0(cancer_type,"_full_cox_table.tsv"), row.names = F, quote = F, sep = "\t")
  
  ### get the transpose
  
  t_keep_matrix = t(keep_matrix)
  t_keep_df = data.frame(keep_locus, t_keep_matrix)
  colnames(t_keep_df) = c("locus", sel_clinical$bcr_patient_barcode)
  
  write.table(t_keep_df, paste0(cancer_type,"_full_cox_table_transpose.tsv"), row.names = F, quote = F, sep = "\t")
  
  return(keep_matrix)
  
}


# setwd("/Users/ginny/Dropbox/PREXA_Rproject/tmp_download/")

#stomach_all = fread("STAD_whole.tsv", stringsAsFactors = F)
#sel_clinical = fread ("stad_clinical.tsv",stringsAsFactors = F)

### focus only on mutations which happen on 10 or more patients

#sel_loci <- stomach_all %>% group_by(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos,mutationaa, aapos, protlen)%>%filter(n()>=10)

###fix some possible bugs
#this_locus_patient <- sel_loci %>% mutate(locus_name = paste(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos, mutationaa,aapos,protlen,sep="-")) %>% ungroup() %>% select(locus_name, sampleid)


#whole_stad_matrix = get_locus_patient_matrix(sel_loci,"stad",sel_clinical)




