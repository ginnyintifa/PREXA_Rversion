### 8 map the domain counts to patient

library(dplyr)
library(magrittr)
library(data.table)


domain_build_tb_each_colum=function(cdf,one_domain)
{
    #one_ptm = ptm_to_map[10765,]
    #cdf = sel_snp
    
    countit=rep(0,nrow(sel_clinical))
    
    this_cdf <- cdf %>% filter (uni_swiss == one_domain$ID | uni_trembl == one_domain$ID)
    this_start_pos <- one_domain$start_position
    this_end_pos <- one_domain$end_position
    
    if (nrow(this_cdf>=0))
    {
        for (i in 1:nrow(this_cdf))
        {
            #i =4943
            cdf_pos=as.numeric(as.character(this_cdf$aapos[i]))
            
            cdf_pos_start=as.numeric(as.character(this_cdf$aaposstart[i]))
            cdf_pos_end=as.numeric(as.character(this_cdf$aaposend[i]))
            
            cdf_patient=this_cdf$sampleid[i]
            
            whichpatient=which(sel_clinical$bcr_patient_barcode==cdf_patient)
            
            if (is.na(cdf_pos)==F &&
            cdf_pos_start>=(this_start_pos)&&
            cdf_pos_end<=(this_end_pos))
            {
                countit[whichpatient]=countit[whichpatient]+1
            }
            
        }
        
    }
    
    countit
    
}



get_patient_domain_mapping = function(sel_clinical,stomach_all, snp_domain_mapping,
output_name, transpose_output_name)
{
    
    
    domain_to_map <- snp_domain_mapping %>% select (protein_ID, ID, start_position, end_position, domain_name, domain_type)
    
    domain_prots_relevant = unique(snp_domain_mapping$ID)
    domain_sel_snp <- stomach_all %>% filter (uni_swiss %in% domain_prots_relevant | uni_trembl %in% domain_prots_relevant)
    
    pos=as.list(domain_sel_snp$aapos)
    
    aaposstart=sapply(pos, function(x) strsplit(x, split="-")[[1]][1])
    aaposend=sapply(pos, function(x) strsplit(x, split="-")[[1]][2])
    
    singlepos=which(is.na(aaposend)==T)
    
    aaposend[singlepos]=aaposstart[singlepos]
    
    new_domain_sel_snp=cbind(domain_sel_snp, aaposstart,aaposend)
    
    
    
    order_patient = sel_clinical$bcr_patient_barcode
    
    
    patient_domain_snp_tb=matrix(0,nrow(sel_clinical),nrow(domain_to_map))
    
    
    for (a in 1:ncol(patient_domain_snp_tb))
    {
        
        patient_domain_snp_tb[,a]=domain_build_tb_each_colum(new_domain_sel_snp,domain_to_map[a,])
        if (a%%1000 ==0)
        {
            cat(a)
            cat("\n")
            
        }
        
    }
    
    
    domain_colnames_matrix=sapply(1:nrow(domain_to_map), function(x) paste(domain_to_map$ID[x],domain_to_map$start_position[x],domain_to_map$end_position[x],domain_to_map$domain_name[x],domain_to_map$domain_type[x],sep="_"))
    
    
    sel_domains=which(apply(patient_domain_snp_tb,2,sum)>0)
    sel_domain_tb=patient_domain_snp_tb[,sel_domains]
    sel_domain_names= domain_colnames_matrix[sel_domains]
    
    sel_domain_df=as.data.frame(sel_domain_tb)
    
    colnames(sel_domain_df)=sel_domain_names
    
    patient_sel_domain_df = data.frame(sel_clinical$bcr_patient_barcode, sel_domain_df)
    
    #output_name = "stad_patient_domain_table.tsv"
    
    write.table(patient_sel_domain_df, output_name, row.names = F, quote = F, sep = "\t")
    
    
    
    
    
    data_domain_trans= t(sel_domain_tb)
    
    data_domain_trans_df= as.data.frame(data_domain_trans)
    
    colnames(data_domain_trans_df)=sel_clinical$bcr_patient_barcode
    
    trans_domain = data.frame(sel_domain_names, data_domain_trans_df)
    write.table(trans_domain, transpose_output_name, row.names = F, quote = F, sep = "\t")
    
    
    
}






#sel_clinical=fread("stad_clinical.tsv",stringsAsFactors = F)
#stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)
#snp_domain_mapping=fread("snp_domain_mapping_sig_stad.tsv", stringsAsFactors = F)

#get_patient_domain_mapping(sel_clinical, stomach_all, snp_domain_mapping,
#"stad_patient_domain_table.tsv",
#"stad_patient_domain_table_transpose.tsv")



