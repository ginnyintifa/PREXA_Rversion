## script to get the data matrix for boxplot


get_df_all = function (silent_coef_p_q,
                       misSense_coef_p_q, 
                       nonSense_coef_p_q, 
                       spliceSite_coef_p_q,
                       inFrameDel_coef_p_q, 
                       inFrameIns_coef_p_q, 
                       frameShiftDel_coef_p_q, 
                       frameShiftIns_coef_p_q,
                       nonStop_coef_p_q,
                       tFlank_coef_p_q,
                       fFlank_coef_p_q,
                       tUTR_coef_p_q,
                       fUTR_coef_p_q,
                       IGR_coef_p_q,
                       intron_coef_p_q,
                       RNA_coef_p_q,
                       targetedRegion_coef_p_q,
                       translationStartSite_coef_p_q,
                       ptm_coef_p_q,
                       domain_coef_p_q)
{
  
  all_snp_data=rbind(silent_coef_p_q,
                     misSense_coef_p_q, 
                     nonSense_coef_p_q, 
                     spliceSite_coef_p_q,
                     inFrameDel_coef_p_q, 
                     inFrameIns_coef_p_q, 
                     frameShiftDel_coef_p_q, 
                     frameShiftIns_coef_p_q,
                     nonStop_coef_p_q,
                     tFlank_coef_p_q,
                     fFlank_coef_p_q,
                     tUTR_coef_p_q,
                     fUTR_coef_p_q,
                     IGR_coef_p_q,
                     intron_coef_p_q,
                     RNA_coef_p_q,
                     targetedRegion_coef_p_q,
                     translationStartSite_coef_p_q
  )
  
  
  all_labels = c(rep("silent",nrow(silent_coef_p_q)),
                 rep("misSense",nrow(misSense_coef_p_q)),
                 rep("nonSense",nrow(nonSense_coef_p_q)),
                 rep("spliceSite",nrow(spliceSite_coef_p_q)),
                 rep("inFrameDel",nrow(inFrameDel_coef_p_q)),
                 rep("inFrameIns",nrow(inFrameIns_coef_p_q)),
                 rep("frameShiftDel",nrow(frameShiftDel_coef_p_q)),
                 rep("frameShiftIns",nrow(frameShiftIns_coef_p_q)),
                 rep("nonStop",nrow(nonStop_coef_p_q)),
                 rep("tFlank",nrow(tFlank_coef_p_q)),
                 rep("fFlank", nrow(fFlank_coef_p_q)),
                 rep("tUTR", nrow(tUTR_coef_p_q)),
                 rep("fUTR", nrow(fUTR_coef_p_q)),
                 rep("IGR", nrow(IGR_coef_p_q)),
                 rep("intron", nrow(intron_coef_p_q)),
                 rep("RNA", nrow(RNA_coef_p_q)),
                 rep("targetedRegion", nrow(targetedRegion_coef_p_q)),
                 rep("translationStartSite",nrow(translationStartSite_coef_p_q)))
  
  
  
  
  
  all_snp_data = data.frame(coef = all_snp_data$coef_value, pValue = all_snp_data$uni_p_value, qValue = p.adjust(all_snp_data$uni_p_value,"BH") ,all_labels)
  
  colnames(all_snp_data) = c("coef", "pValue","qValue","all_labels")
  
  
  all_pro_data = rbind (ptm_coef_p_q, domain_coef_p_q)
  all_labels = c(rep("ptm", nrow(ptm_sel_coef_p_q)),
                 rep("domain", nrow(domain_sel_coef_p_q)))
  
  all_pro_data = data.frame(coef = all_pro_data$coef_value, pValue = all_pro_data$uni_p_value, qValue = all_pro_data$q_value ,all_labels)
  
  colnames(all_pro_data) = c("coef", "pValue","qValue","all_labels")
  
  all_data = rbind(all_snp_data, all_pro_data)
  
  ###get the data processed so that data with less than 100 row are not there
  
  
  sel_all_data <- all_data %>% 
    group_by(all_labels) %>% 
    mutate (the_count = n() ) %>%
    filter (the_count >= 100)
  
  #gs <- sel_all_data %>% group_by(all_labels) %>% summarise(n())
  
  
  arranged_names <- sel_all_data %>% 
    group_by(all_labels) %>% 
    summarize(median_q = median(qValue)) %>%
    arrange (median_q) %>%
    select(all_labels)
  
  
  sel_all_data$all_labels <- factor(sel_all_data$all_labels, levels = arranged_names$all_labels )
  
  
  write.table(sel_all_data, "sel_all_coef_p_q.tsv", quote= F, row.names = F,
              sep = "\t")
  
   
  return (sel_all_data)
  
  
}

