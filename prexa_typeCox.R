###  11 Get cox model built on ptm/domain/prot
library(survival)
library(dplyr)
library(magrittr)
library(data.table)

### I should be able to know which row is which prot/ptm/domain 
### now I don't know...




get_type_q_value = function(patient_type_data, data_type)
{
    
    data_outof_table = as.matrix(patient_type_data[,-1])
    
    
    sum_type_counts = apply(data_outof_table,2,sum)
    
    
    which_sel = which(sum_type_counts>=10)
    
    sel_type = data_outof_table[,which_sel]
    
    sel_df = data.frame(patient_id = patient_type_data[,1],sel_type)
    
    colnames(sel_df)[1]="patient_id"
    
    write.table(sel_df, paste(data_type,"cox_table.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
    
    unit_identifier = colnames(sel_type)
    t_sel_type = t(sel_type)
    
    t_sel_df = data.frame(units = unit_identifier,t_sel_type)
    
    colnames(t_sel_df) = c("units", patient_type_data[,1])
    
    write.table(t_sel_df, paste(data_type,"cox_table_transpose.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
    
    
    ### get the identifier
    
    ### up to here is correct because the clinical info is not linked yet
    
    
    if (ncol(sel_type)>0)
    {
        
        identifier = colnames(sel_type)
        
        uni_p_value = rep(0, ncol(sel_type))
        coef_value = rep(0, ncol(sel_type))
        
        
        ### need to sort the rows of sel_type in the order of sorted_clinical info
        
        
        
        
        for (p in 1: ncol(sel_type))
        {
            
            
            this_cox = coxph(Surv(time, status) ~ sel_type[,p] +
            gender + race + age )
            sumit= summary(this_cox)
            
            uni_p_value[p] = sumit$coefficients[1,5]
            coef_value[p] = sumit$coefficients[1,1]
            
        }
        
        q_value = p.adjust(uni_p_value, method = "BH")
        
        
        #### try to get the identifier of each row 
        
        cox_out = cbind(identifier, coef_value, uni_p_value, q_value)
        
        write.table(cox_out, paste(data_type,"cox_out.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
        
        
        
        
        
        return(cox_out)
    }else{
        return (matrix(0,0,3))
    }
    
    
    
}


