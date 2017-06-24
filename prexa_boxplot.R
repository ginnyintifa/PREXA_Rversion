### 13 Get boxplot with ggplot2
library(ggplot2)
library(dplyr)
library(magrittr)
library(data.table)




get_meta_boxplot = function(silent_coef_p_q,
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
)
{
    
    all_data=rbind(silent_coef_p_q,
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
    rep("translationStartSite", nrow(translationStartSite_coef_p_q)))
    
    
    
    
    
    all_data = data.frame(coef = all_data[,1], pValue = all_data[,2], qValue = all_data[,3],all_labels)
    
    
    arranged_names_pValue <- all_data %>%
      group_by(all_labels) %>%
      summarize(median_p = median(pValue)) %>%
      arrange (median_p) %>%
      select(all_labels)
    
    p_all_data = all_data
    
    p_all_data$all_labels <- factor(p_all_data$all_labels, levels = arranged_names_pValue$all_labels )
    
    ggplot(data = p_all_data, aes(x = all_labels, y = -log10(pValue)))+
    geom_boxplot()+
    ggtitle ("boxplot of -log10(pValue)")
    
    
    ###
    arranged_names_qValue <- all_data %>%
      group_by(all_labels) %>%
      summarize(median_q = median(qValue)) %>%
      arrange (median_q) %>%
      select(all_labels)
    
    q_all_data = all_data
    
    q_all_data$all_labels <- factor(q_all_data$all_labels, levels = arranged_names_qValue$all_labels )
    


    ggplot(data = q_all_data, aes(x = all_labels, y = -log10(qValue)))+
    geom_boxplot()+
    ggtitle ("boxplot of -log10(qValue)")
    
    
    
    ###
    
    arranged_names_coefValue <- all_data %>%
      group_by(all_labels) %>%
      summarize(median_coef = median(coef)) %>%
      arrange (desc(median_coef)) %>%
      select(all_labels)
    
    coef_all_data = all_data
    
    coef_all_data$all_labels <- factor(coef_all_data$all_labels, levels = arranged_names_coefValue$all_labels )
    
    
    
    ggplot(data = coef_all_data, aes(x = all_labels, y = exp(coef)))+
    geom_boxplot()+
    ggtitle ("boxplot of exp(coefficient)")
    
    

    
}


