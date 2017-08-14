
library(survival)

get_this_column = function(i, all_locus, this_locus_patient, sel_clinical)
{
    
    this_column = rep(0, length(sel_clinical$bcr_patient_barcode))
    get_patient <- this_locus_patient %>% filter(locus_name == all_locus[i])%>% select(sampleid)
    find_patient = which(sel_clinical$bcr_patient_barcode %in% get_patient$sampleid)
    this_column[find_patient]=1
    
    return (this_column)
    
}


get_locus_patient_df = function (this_type, sel_loci, sel_clinical)
{
    
    this_data <- sel_loci%>%filter(mutationtype == this_type)
    
    locus_name =  paste(this_data$genename, this_data$uni_swiss, this_data$uni_trembl,this_data$chromosome, this_data$mutationtype, this_data$mutationpos,this_data$mutationaa,this_data$aapos,
    this_data$protlen,sep="-")
    
    #this_locus_patient <- this_data %>% mutate(locus_name = paste(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos, mutationaa,aapos,protlen,sep="-")) %>% ungroup() %>% select(locus_name, sampleid)
    this_locus_patient = data.frame(locus_name, sampleid = this_data$sampleid)
    
    #this_matrix = matrix(0, length(right_patient_info$bcr_patient_barcode), length(all_locus))
    
    all_locus = unique(locus_name)
    the_matrix = sapply(1:length(all_locus), function(i){
        it = get_this_column(i, all_locus, this_locus_patient,sel_clinical )
        return(it)
    } )
    
    the_matrix_keep_mutation = which(apply(the_matrix,2,sum)>=10)
    
    keep_matrix = the_matrix[,the_matrix_keep_mutation]
    keep_locus = all_locus[the_matrix_keep_mutation]
    
    
    keep_df = data.frame(sel_clinical$bcr_patient_barcode,keep_matrix)
    colnames(keep_df) = c("patient_id",keep_locus)
    
    write.table(keep_df, paste(this_type,"cox_table.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
    
    ### so this can be passed to the next step
    
    return(keep_df)
    
}

### get each locus to cox model and get the p/q value for them

get_q_value = function (this_type, sel_loci, sel_clinical)

{
    
    
    sel_get = get_locus_patient_df(this_type, sel_loci, sel_clinical)
    
    identifier = colnames(sel_get)[-1]
    sel_table = as.matrix(sel_get[,-1])
    
    
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
        
        ### trying to add identifier
        cox_out = cbind(identifier, coef_value, uni_p_value, q_value)
        
        write.table(cox_out, paste(this_type,"cox_out.tsv",sep="_"), row.names = F, quote = F, sep = "\t")
        
        
        return(cox_out)
        
    }else(
    return (matrix(0,0,3))
    )
    
    
    
}




