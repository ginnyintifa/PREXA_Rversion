
### 9 Get the count of each proteins, regardless of PTM or domain
library(dplyr)
library(magrittr)
library(data.table)

protein_build_tb_each_colum=function(cdf,one_protein)
{
    countit=rep(0,nrow(sel_clinical))
    
    this_cdf <- cdf %>% filter (uni_swiss == one_protein | uni_trembl == one_protein)
    
    if (nrow(this_cdf>=0))
    {
        for (i in 1:nrow(this_cdf))
        {
            
            
            cdf_patient=this_cdf$sampleid[i]
            
            whichpatient=which(sel_clinical$bcr_patient_barcode==cdf_patient)
            
            
            countit[whichpatient]=countit[whichpatient]+1
            
            
        }
        
    }
    
    countit
    
}




get_patient_prot_mapping = function (sel_clinical, stomach_all,
output_name, transpose_output_name)
{
    sel_proteins <- stomach_all%>%select(uni_swiss, uni_trembl)
    
    unique_proteins = unique( c(sel_proteins$uni_swiss, sel_proteins$uni_trembl))
    ###need to get rid of the blank proteins
    
    uni_prots = unique_proteins[which(nchar(unique_proteins)>0)]
    
    order_patient = sel_clinical$bcr_patient_barcode
    
    
    #cdf = tb_stomach_all
    #one_protein = sel_protein[1]
    
    patient_protein_snp_tb=matrix(0,nrow(sel_clinical),length(uni_prots))
    
    
    for (a in 1:ncol(patient_protein_snp_tb))
    {
        
        patient_protein_snp_tb[,a]=protein_build_tb_each_colum(stomach_all,
        uni_prots[a])
        if (a%%1000 ==0)
        {
            cat(a)
            cat("\n")
            
        }
        
    }
    
    
    
    
    sel_proteins_which=which(apply(patient_protein_snp_tb,2,sum)>0)
    sel_proteins_tb=patient_protein_snp_tb[,sel_proteins_which]
    sel_proteins_names= uni_prots[sel_proteins_which]
    
    sel_proteins_df=as.data.frame(sel_proteins_tb)
    
    colnames(sel_proteins_df)=sel_proteins_names
    
    patient_sel_proteins_df = data.frame(sel_clinical$bcr_patient_barcode, sel_proteins_df)
    
    write.table(patient_sel_proteins_df, output_name, row.names = F, quote = F, sep = "\t")
    
    
    
    
    
    data_proteins_trans= t(sel_proteins_tb)
    
    data_proteins_trans_df= as.data.frame(data_proteins_trans)
    
    colnames(data_proteins_trans_df)=sel_clinical$bcr_patient_barcode
    
    trans_proteins = data.frame(sel_proteins_names, data_proteins_trans_df)
    write.table(trans_proteins, transpose_output_name, row.names = F, quote = F, sep = "\t")
    
}



#sel_clinical=fread("stad_clinical.tsv",stringsAsFactors = F)
#stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)

#get_patient_prot_mapping(sel_clinical, stomach_all,
#"stad_patient_prot_table.tsv",
#"stad_patient_prot_table_transpose.tsv")
