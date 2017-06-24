

# 10 SNP mutation type model building

library(survival)
library(dplyr)
library(magrittr)
library(data.table)


get_this_column = function(i, all_locus, this_locus_patient)
{
    
    
    this_column = rep(0, length(sel_clinical$bcr_patient_barcode))
    get_patient <- this_locus_patient %>% filter(locus_name == all_locus[i])%>% select(sampleid)
    find_patient = which(sel_clinical$bcr_patient_barcode %in% get_patient$sampleid)
    this_column[find_patient]=1
    
    return (this_column)
    
}


get_locus_patient_matrix = function (this_type, sel_loci)
{
    
    this_data <- sel_loci%>%filter(mutationtype == this_type)
    
    this_locus_patient <- this_data %>% mutate(locus_name = paste(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos, mutationaa,aapos,protlen,sep="-")) %>% ungroup() %>% select(locus_name, sampleid)
    
    
    #this_matrix = matrix(0, length(right_patient_info$bcr_patient_barcode), length(all_locus))
    
    all_locus = unique(this_locus_patient$locus_name)
    the_matrix = sapply(1:length(all_locus), function(i){
        it = get_this_column(i, all_locus, this_locus_patient)
        return(it)
    } )
    
    the_matrix_keep_mutation = which(apply(the_matrix,2,sum)>=10)
    
    keep_matrix = the_matrix[,the_matrix_keep_mutation]
    keep_locus = all_locus[the_matrix_keep_mutation]
    
    
    keep_df = data.frame(sel_clinical$bcr_patient_barcode,keep_matrix)
    colnames(keep_df) = c("patient_id",keep_locus)
    
    write.table(keep_df, paste(this_type,"cox_table.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
    
    return(keep_matrix)
    
}


get_q_value = function (this_type, sel_loci)

{
    
    # this_type = "Targeted_Region"
    #sel_loci <- stomach_all %>% group_by(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos,mutationaa, aapos, protlen)%>%filter(n()>=10)
    
    sel_table = as.matrix(get_locus_patient_matrix(this_type, sel_loci))
    
    #sel_ptm = silent_cox_table
    uni_p_value = rep(0, ncol(sel_table))
    coef_value = rep(0, ncol(sel_table))
    
    if (ncol(sel_table)>0)
    {
        for (p in 1: ncol(sel_table))
        {
            
            
            this_cox = coxph(Surv(time, status) ~ sel_table[,p] +
            gender + race + age )
            sumit= summary(this_cox)
            
            uni_p_value[p] = sumit$coefficients[1,5]
            coef_value[p] = sumit$coefficients[1,1]
            
        }
        
        q_value = p.adjust(uni_p_value, method = "BH")
        
        
        cox_out = cbind(coef_value, uni_p_value, q_value)
        
        write.table(cox_out, paste(this_type,"cox_out.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
        
        
        return(cox_out)
        
    }else(
    return (matrix(0,0,3))
    )
    
    
    
}


#setwd("/Users/ginny/tcga_june_gastric_cancer/TCGA_STAD/")
#stomach_all = fread("STAD_whole.tsv", stringsAsFactors = F)
#sel_clinical = fread("stad_clinical.tsv", stringsAsFactors = F)

    gender = sel_clinical$gender
    race = sel_clinical$race
    age = sel_clinical$age_at_initial_pathologic_diagnosis
    
    not_av_age = which(substr(age, 1,1)=="[")
    age[not_av_age] = mean(as.numeric(age[-not_av_age]))
    age = as.numeric(age)
    
    not_av_race = which(substr(race, 1,1)=="[")
    race[not_av_race] = "OTHER"
    race = as.factor(race)
    
    gender = as.factor(gender)
    
    
    time = sel_clinical$survival_time
    status = sel_clinical$survival_ind




sel_loci <- stomach_all %>% group_by(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos,mutationaa, aapos, protlen)%>%filter(n()>=10)


silent_coef_p_q = get_q_value("Silent", sel_loci)
misSense_coef_p_q = get_q_value("Missense_Mutation",sel_loci)
nonSense_coef_p_q = get_q_value("Nonsense_Mutation", sel_loci)
spliceSite_coef_p_q= get_q_value("Splice_Site", sel_loci)
inFrameDel_coef_p_q = get_q_value("In_Frame_Del", sel_loci)
inFrameIns_coef_p_q = get_q_value("In_Frame_Ins", sel_loci)
frameShiftDel_coef_p_q = get_q_value("Frame_Shift_Del", sel_loci)
frameShiftIns_coef_p_q = get_q_value("Frame_Shift_Ins", sel_loci)
nonStop_coef_p_q = get_q_value("Nonstop_Mutation", sel_loci)
tFlank_coef_p_q = get_q_value("3'Flank", sel_loci)
fFlank_coef_p_q = get_q_value("5'Flank", sel_loci)
tUTR_coef_p_q = get_q_value("3'UTR", sel_loci)
fUTR_coef_p_q = get_q_value("5'UTR", sel_loci)
IGR_coef_p_q = get_q_value("IGR", sel_loci)
intron_coef_p_q = get_q_value("Intron", sel_loci)
RNA_coef_p_q = get_q_value("RNA", sel_loci)
targetedRegion_coef_p_q = get_q_value("Targeted_Region", sel_loci)
translationStartSite_coef_p_q = get_q_value("Translation_Start_Site", sel_loci)
