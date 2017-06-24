###  11 Get cox model built on ptm/domain/prot
library(survival)
library(dplyr)
library(magrittr)
library(data.table)


get_type_q_value = function(patient_type_data, data_type)
{
    data_outof_table = as.matrix(patient_type_data[,-1])
    
    sum_type_counts = apply(data_outof_table,2,sum)
    
    
    which_sel = which(sum_type_counts>=10)
    
    sel_type = data_outof_table[,which_sel]
    
    
    write.table(sel_type, paste(data_type,"cox_table.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
    
    
    sel_type_id = colnames(sel_type)
    
    if (ncol(sel_type)>0)
    {
        uni_p_value = rep(0, ncol(sel_type))
        coef_value = rep(0, ncol(sel_type))
        
        
        for (p in 1: ncol(sel_type))
        {
            
            
            this_cox = coxph(Surv(time, status) ~ sel_type[,p] +
            gender + race + age )
            sumit= summary(this_cox)
            
            uni_p_value[p] = sumit$coefficients[1,5]
            coef_value[p] = sumit$coefficients[1,1]
            
        }
        
        q_value = p.adjust(uni_p_value, method = "BH")
        
        cox_out = cbind(coef_value, uni_p_value, q_value)
        
        write.table(cox_out, paste(data_type,"cox_out.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
        
        
        return(cox_out)
    }else{
        return (matrix(0,0,3))
    }
    
}


#STAD

stad_patient_ptm = fread("stad_patient_ptm_table.tsv", stringsAsFactors = F, header = T)

stad_patient_domain =fread("stad_patient_domain_table.tsv", stringsAsFactors = F, header = T)

stad_patient_prot =fread("stad_patient_prot_table.tsv", stringsAsFactors = F, header = T)

sel_clinical = fread("stad_clinical.tsv",stringsAsFactors = F)


gender =sel_clinical$gender
race = sel_clinical$race
age = sel_clinical$age_at_initial_pathologic_diagnosis

not_av_age = which(substr(age, 1,1)=="[")
age[not_av_age] = mean(as.numeric(age[-not_av_age]))
age = as.numeric(age)

not_av_race = which(substr(race, 1,1)=="[")
race[not_av_race] = "OTHER"
race = as.factor(race)

gender = as.factor(gender)



ptm_coef_p_q = get_type_q_value(stad_patient_ptm,"ptm")
domain_coef_p_q = get_type_q_value(stad_patient_domain,"domain")
prot_coef_p_q = get_type_q_value(stad_patient_prot,"prot")

