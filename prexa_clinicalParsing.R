


## 6 Parse the new patient information


library(dplyr)
library(magrittr)
library(data.table)



get_clinical = function (tcga_follow_up, tcga_patient, output_name)
{
    
    tcga_follow_info <- tcga_follow_up %>%
    select(bcr_patient_barcode,vital_status,
    last_contact_days_to, death_days_to,
    new_tumor_event_dx_indicator,new_tumor_event_dx_indicator)%>%
    filter(vital_status=="Alive" | vital_status =="Dead")
    ###get rid of the not availabel cases already
    
    ## get all patients
    my_patients = unique(tcga_follow_info$bcr_patient_barcode)
    
    p_ind = rep(0,length(my_patients))
    
    p_time = rep(0, length(my_patients))
    
    
    for (i in 1:length(my_patients))
    {
        #i=2
        
        cat(i)
        cat("\n")
        
        which_row = tcga_follow_info$bcr_patient_barcode == my_patients[i]
        
        all_ind = tcga_follow_info$vital_status[which_row]
        all_time = c(tcga_follow_info$last_contact_days_to[which_row],
        tcga_follow_info$death_days_to[which_row])
        
        
        if ("Dead" %in% all_ind)
        {
            p_ind[i] = 1
        }
        
        if(all(substr(all_time,1,1)=="["))
        {
            p_time[i] = NA
        }else
        {
            numeric_all_time = as.numeric(all_time[which(substr(all_time,1,1) != "[")])
            
            p_time[i] = max(numeric_all_time)
            
        }
        
        
        
    }
    
    
    which_ava_patients = which(is.na(p_time) ==F & p_time >0)
    
    
    ava_patients = my_patients[which_ava_patients]
    ava_p_ind = p_ind[which_ava_patients]
    ava_p_time = p_time[which_ava_patients]
    
    patients_survival = data.frame(ava_patients, ava_p_ind, ava_p_time)
    
    
    
    
    
    tcga_patient_info <- tcga_patient %>%
    select(bcr_patient_barcode, tumor_grade, gender, race, ethnicity,
    age_at_initial_pathologic_diagnosis, anatomic_neoplasm_subdivision,
    tumor_grade, ajcc_pathologic_tumor_stage,tumor_status,
    family_history_of_stomach_cancer, number_of_relatives_with_stomach_cancer,
    targeted_molecular_therapy,radiation_treatment_adjuvant,
    tumor_sample_procurement_country,initial_pathologic_dx_year,
    lymph_nodes_examined, residual_tumor
    ) %>%
    filter(bcr_patient_barcode %in% ava_patients)
    
    as_patient_order <- patients_survival %>%
    arrange(tcga_patient_info$bcr_patient_barcode)
    
    tcga_clinical <- tcga_patient_info %>%
    mutate(survival_ind = as_patient_order$'ava_p_ind',
    survival_time = as_patient_order$'ava_p_time')%>%
    arrange(survival_ind)
    
    # output_name = "tcga_stad_clinical.tsv"
    
    write.table(tcga_clinical, output_name, sep= "\t", row.names = F, quote = F)
    
    
}


### for stomach cancer

#setwd("/Users/ginny/tcga_june_gastric_cancer/TCGA_STAD")

#tcga_follow_up <- fread("tcga_follow_up.txt", stringsAsFactors = F)
#tcga_patient <- fread("tcga_patient.txt", stringsAsFactors = F)
#tcga_stad_patient_info = get_clinical(tcga_follow_up, tcga_patient, "stad_clinical.tsv")


