
## 5 Map the variance to domains, and get the counts

### consider only mappable mutations





#domain_data=fread("domain_sorted_df.tsv", header = T, stringsAsFactors = F)

library(dplyr)
library(magrittr)
library(data.table)



treat_each_row_domain=function(proteinID, start_pos,end_pos, whole_snp)
{
    a_total_mut=0
    a_missense_mut=0
    a_nonsense_mut=0
    a_silent_mut=0
    a_splice_mut=0
    a_frame_shift_ins=0
    a_frame_shift_del=0
    a_in_frame_ins=0
    a_in_frame_del=0
    
    
    wanted_snp_row=which(whole_snp$uni_swiss==proteinID|whole_snp$uni_trembl==proteinID)
    snp_wanted=whole_snp[wanted_snp_row,]
    if (length(wanted_snp_row)!=0)
    {
        
        near_mut_start=which(as.numeric(as.character(snp_wanted$aaposstart))<=(end_pos)&as.numeric(as.character(snp_wanted$aaposstart))>=(start_pos))
        near_mut_end=which(as.numeric(as.character(snp_wanted$aaposend))<=(end_pos)&as.numeric(as.character(snp_wanted$aaposend))>=(start_pos))
        
        near_mut_which=unique(c(near_mut_start,near_mut_end))
        a_total_mut=a_total_mut+length(near_mut_which)
        
        snp_mut_near=snp_wanted[near_mut_which,]
        
        missense_which=which(snp_mut_near$mutationtype=="Missense_Mutation")
        nonsense_which=which(snp_mut_near$mutationtype=="Nonsense_Mutation")
        silent_which=which(snp_mut_near$mutationtype=="Silent")
        splice_which=which(snp_mut_near$mutationtype=="Splice_Site")
        frame_shift_ins_which=which(snp_mut_near$mutationtype=="Frame_Shift_Ins")
        frame_shift_del_which=which(snp_mut_near$mutationtype=="Frame_Shift_Del")
        in_frame_ins_which=which(snp_mut_near$mutationtype=="In_Frame_Ins")
        in_frame_del_which=which(snp_mut_near$mutationtype=="In_Frame_Del")
        
        
        
        a_missense_mut=a_missense_mut+length(missense_which)
        a_nonsense_mut=a_nonsense_mut+length(nonsense_which)
        a_silent_mut=a_silent_mut+length(silent_which)
        a_splice_mut=a_splice_mut+length(splice_which)
        a_frame_shift_ins=a_frame_shift_ins+length(frame_shift_ins_which)
        a_frame_shift_del=a_frame_shift_del+length(frame_shift_del_which)
        a_in_frame_ins=a_in_frame_ins+length(in_frame_ins_which)
        a_in_frame_del=a_in_frame_del+length(in_frame_del_which)
        
    }
    
    
    all_mut_info=c( a_total_mut,a_missense_mut,a_nonsense_mut,a_silent_mut, a_splice_mut,
    a_frame_shift_ins,a_frame_shift_del,a_in_frame_ins,a_in_frame_del)
    
    
    all_mut_info
    
    
}





get_domain_mapping = function(stomach_all, output_name)
{
    pos=as.list(stomach_all$aapos)
    
    aaposstart=sapply(pos, function(x) strsplit(x, split="-")[[1]][1])
    aaposend=sapply(pos, function(x) strsplit(x, split="-")[[1]][2])
    
    singlepos=which(is.na(aaposend)==T)
    
    aaposend[singlepos]=aaposstart[singlepos]
    
    new_snp_df=cbind(stomach_all, aaposstart,aaposend)
    
    
    cancer_snp_domain=sapply(1:nrow(domain_data),
    function(i){
        protid=domain_data$ID[i]
        start_pos=as.numeric(domain_data$start_position[i])
        end_pos=as.numeric(domain_data$end_position[i])
        
        this_mut=treat_each_row_domain(protid, start_pos, end_pos, new_snp_df)
        
        
        if(i%%10000==0)
        {
            cat(i)
            cat("\n")
        }
        
        return(this_mut)
        
        
        
    })
    
    
    t_cancer_snp_domain=t(cancer_snp_domain)
    
    
    
    colnames(t_cancer_snp_domain)=c( "total","missense","nonsense","silent","splice_site","frame_shift_ins","frame_shift_del","in_frame_ins","in_frame_del")
    
    t_cancer_snp_domain_sig=t_cancer_snp_domain[which(apply(t_cancer_snp_domain,1,sum)>0),]
    
    domain_data_sig=domain_data[which(apply(t_cancer_snp_domain,1,sum)>0),]
    
    #output_name = "snp_domain_mapping_sig.tsv"
    
    write.table(data.frame(domain_data_sig,t_cancer_snp_domain_sig),output_name ,row.names = F,quote=F,sep="\t")
    
}




#domain_data=fread("/Users/ginny/tcga_june_gastric_cancer/domain_sorted_df.tsv", header = T, stringsAsFactors = F)


####1:42
#stomach_all=fread("STAD_whole.tsv", header = T, stringsAsFactors = F)
#get_domain_mapping (stomach_all, "snp_domain_mapping_sig.tsv")

