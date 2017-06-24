
## 7 map the PTM counts to patient
library(dplyr)
library(magrittr)
library(data.table)




build_tb_each_colum=function(cdf,one_ptm)
{
    #one_ptm = ptm_to_map[10765,]
    #cdf = sel_snp
    
    countit=rep(0,nrow(sel_clinical))
    
    this_cdf <- cdf %>% filter (uni_swiss == one_ptm$ID | uni_trembl == one_ptm$ID)
    this_pos <- one_ptm$position
    
    if (nrow(this_cdf>=0))
    {
        for (i in 1:nrow(this_cdf))
        {
            
            cdf_pos = as.numeric(as.character(this_cdf$aapos[i]))
            cdf_pos_start=as.numeric(as.character(this_cdf$aaposstart[i])) ##this line may be problemmatic when it is an area(with "-"), I need to create the same start and end way to make it correct.do this before hand
            cdf_pos_end=as.numeric(as.character(this_cdf$aaposend[i]))
            
            cdf_patient=this_cdf$sampleid[i]
            
            whichpatient=which(sel_clinical$bcr_patient_barcode==cdf_patient)
            
            if (is.na(cdf_pos)==F && cdf_pos_start>=(this_pos-5)&&cdf_pos_end<=(this_pos+5))
            {
                countit[whichpatient]=countit[whichpatient]+1
            }
            
        }
        
    }
    
    countit
    
}



get_patient_ptm_mapping = function (sel_clinical, stomach_all, snp_ptm_mapping,
output_name, transpose_output_name)
{
    
    ptm_to_map <- snp_ptm_mapping %>% select (protein_ID, ID, position, ptm_type)
    
    prots_relevant = unique(snp_ptm_mapping$ID)
    sel_snp <- stomach_all %>% filter (uni_swiss %in% prots_relevant | uni_trembl %in% prots_relevant)
    
    pos=as.list(sel_snp$aapos)
    
    aaposstart=sapply(pos, function(x) strsplit(x, split="-")[[1]][1])
    aaposend=sapply(pos, function(x) strsplit(x, split="-")[[1]][2])
    
    singlepos=which(is.na(aaposend)==T)
    
    aaposend[singlepos]=aaposstart[singlepos]
    
    new_sel_snp=cbind(sel_snp, aaposstart,aaposend)
    
    
    
    
    #order_patient = sel_clinical$bcr_patient_barcode
    
    ###order mypatients in the order of live and deceased ?
    
    patient_ptm_snp_tb=matrix(0,nrow(sel_clinical),nrow(ptm_to_map))
    
    for (a in 1:ncol(patient_ptm_snp_tb))
    {
        
        ### I need to check it there is any bug here.
        ### because now sel_snp may contain many entries which have no
        ### position info
        
        ### I need to have the same mapping scheme as the previous ptm-protein mapping
        
        
        patient_ptm_snp_tb[,a]=build_tb_each_colum(new_sel_snp,ptm_to_map[a,])
        
        
        if(a%%1000==0)
        {
            cat(a)
            cat("\n")
        }
        
    }
    
    
    ###make the column names for the matrix, it should be the ptm sites
    
    colnames_matrix=sapply(1:nrow(ptm_to_map), function(x) paste(ptm_to_map$ID[x],ptm_to_map$position[x],ptm_to_map$ptm_type[x],sep="_"))
    
    
    sel_ptms=which(apply(patient_ptm_snp_tb,2,sum)>0)
    sel_ptm_tb=patient_ptm_snp_tb[,sel_ptms]
    sel_ptm_names= colnames_matrix[sel_ptms]
    
    sel_ptm_df=as.data.frame(sel_ptm_tb)
    
    colnames(sel_ptm_df)=sel_ptm_names
    
    patient_sel_ptm_df = data.frame(sel_clinical$bcr_patient_barcode, sel_ptm_df)
    
    #output_name= "stad_patient_ptm_table.tsv"
    
    write.table(patient_sel_ptm_df, output_name, row.names = F, quote = F, sep = "\t")
    
    
    
    
    data_ptm_trans= t(sel_ptm_tb)
    
    data_ptm_trans_df= as.data.frame(data_ptm_trans)
    
    colnames(data_ptm_trans_df)=sel_clinical$bcr_patient_barcode
    
    trans_ptm = data.frame(sel_ptm_names, data_ptm_trans_df)
    
    #transpose_output_name = "stad_patient_ptm_table_transpose.tsv"
    
    write.table(trans_ptm, transpose_output_name, row.names = F, quote = F, sep = "\t")
    
    
}

#sel_clinical=fread("stad_clinical.tsv",stringsAsFactors = F)
#stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)
#snp_ptm_mapping=fread("snp_ptm_mapping_sig_stad.tsv", stringsAsFactors = F)

#get_patient_ptm_mapping(sel_clinical, stomach_all, snp_ptm_mapping,
#"stad_patient_ptm_table.tsv",
#"stad_patient_ptm_table_transpose.tsv")


