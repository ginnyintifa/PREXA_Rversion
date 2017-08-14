### visualize with heatmap on significant genes
### figure out loci and types which are significat after univariate cox regression
### and get heatmap on those


###  I need to figure out the significant ones first from my boxplot data and 
### find the data in the patient mapping matrix and make the heatmap


### ok now get the mutations which are significant


silent_coef_p_q = fread("Silent_cox_out.tsv")
misSense_coef_p_q = fread("Missense_Mutation_cox_out.tsv")
nonSense_coef_p_q = fread("Nonsense_Mutation_cox_out.tsv")
spliceSite_coef_p_q = fread("Splice_Site_cox_out.tsv")
inFrameDel_coef_p_q = fread ("In_Frame_Del_cox_out.tsv")
inFrameIns_coef_p_q = fread ("In_Frame_Ins_cox_out.tsv")
frameShiftDel_coef_p_q = fread ("Frame_Shift_Del_cox_out.tsv")
frameShiftIns_coef_p_q = fread ("Frame_Shift_Ins_cox_out.tsv")
nonStop_coef_p_q = fread ("Nonstop_Mutation_cox_out.tsv")
tFlank_coef_p_q = fread ("3'Flank_cox_out.tsv")
fFlank_coef_p_q = fread ("5'Flank_cox_out.tsv")
tUTR_coef_p_q = fread ("3'UTR_cox_out.tsv")
fUTR_coef_p_q = fread ("5'UTR_cox_out.tsv")
IGR_coef_p_q = fread ("IGR_cox_out.tsv")
intron_coef_p_q = fread ("Intron_cox_out.tsv")
RNA_coef_p_q = fread ("RNA_cox_out.tsv")
targetedRegion_coef_p_q = fread ("Targeted_Region_cox_out.tsv")
translationStartSite_coef_p_q = fread ("Translation_Start_Site_cox_out.tsv")
ptm_coef_p_q = fread ("ptm_cox_out.tsv")
domain_coef_p_q = fread ("domain_cox_out.tsv")
prot_coef_p_q = fread ("prot_cox_out.tsv")


silent_cox = fread("Silent_cox_table.tsv")
misSense_cox = fread("Missense_Mutation_cox_table.tsv")
nonSense_cox = fread("Nonsense_Mutation_cox_table.tsv")
spliceSite_cox = fread("Splice_Site_cox_table.tsv")
inFrameDel_cox = fread ("In_Frame_Del_cox_table.tsv")
inFrameIns_cox= fread ("In_Frame_Ins_cox_table.tsv")
frameShiftDel_cox = fread ("Frame_Shift_Del_cox_table.tsv")
frameShiftIns_cox = fread ("Frame_Shift_Ins_cox_table.tsv")
nonStop_cox = fread ("Nonstop_Mutation_cox_table.tsv")
tFlank_cox = fread ("3'Flank_cox_table.tsv")
fFlank_cox = fread ("5'Flank_cox_table.tsv")
tUTR_cox = fread ("3'UTR_cox_table.tsv")
fUTR_cox = fread ("5'UTR_cox_table.tsv")
IGR_cox = fread ("IGR_cox_table.tsv")
intron_cox = fread ("Intron_cox_table.tsv")
RNA_cox = fread ("RNA_cox_table.tsv")
targetedRegion_cox = fread ("Targeted_Region_cox_table.tsv")
translationStartSite_cox = fread ("Translation_Start_Site_cox_table.tsv")
ptm_cox = fread ("ptm_cox_table.tsv")
domain_cox = fread ("domain_cox_table.tsv")
prot_cox = fread ("prot_cox_table.tsv")






get_data_for_heatmap = function(silent_coef_p_q,
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
                                silent_cox, 
                                misSense_cox, 
                                nonSense_cox, 
                                spliceSite_cox,
                                inFrameDel_cox,
                                inFrameIns_cox,
                                frameShiftDel_cox,
                                frameShiftIns_cox,
                                nonStop_cox,
                                tFlank_cox,
                                fFlank_cox, 
                                tUTR_cox,
                                fUTR_cox,
                                IGR_cox,
                                intron_cox,
                                targetedRegion_cox,
                                translationStartSite_cox
                                )
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
                     translationStartSite_coef_p_q)
  
  all_snp_data$q_value = p.adjust(all_snp_data$uni_p_value)
  
  all_snp_sig <- all_snp_data %>% filter (q_value <= 0.2) %>% select(identifier)
  
  ptm_sig <- ptm_coef_p_q %>% filter (q_value <= 0.2) %>% select(identifier)
  domain_sig <- domain_coef_p_q %>% filter (q_value <= 0.2) %>% select(identifier)
  prot_sig <- prot_coef_p_q %>% filter (q_value <= 0.2) %>% select(identifier)
  
  
  tog = cbind(silent_cox[,-1], misSense_cox[,-1], nonSense_cox[,-1], 
              spliceSite_cox[,-1], inFrameDel_cox[,-1], 
              inFrameIns_cox[,-1], frameShiftDel_cox[,-1], frameShiftIns_cox[,-1],
              nonStop_cox[,-1],tFlank_cox[,-1], fFlank_cox[,-1], 
              tUTR_cox[,-1], fUTR_cox[,-1], IGR_cox[,-1], intron_cox[,-1],
              targetedRegion_cox[,-1], translationStartSite_cox[,-1])
  
  to_sig = tog %>% select (one_of(all_snp_sig$identifier))
  
  trans = t(to_sig) 
  
  colnames(trans) = silent_cox$patient_id
  
  stad_sel_sort_clinical = fread("stad_sel_sort_clinical.tsv",stringsAsFactors = F)
  
  sel_patient = stad_sel_sort_clinical$bcr_patient_barcode
  
  trans_order = trans[,sel_patient]
  
  write.table(trans_order, "significant_snp.tsv",quote= F, row.names = F,
              sep = "\t")
  
  #library(gplots)
  
  
  #cb = c(-0.1,0,1.1)
  
  #heatmap.2(trans_order, trace = "n",Rowv = T, Colv= F,
  #main ="snp", dendrogram = "row", col=colorRampPalette(c("white","black")),breaks=cb, margins = c(2,14))
  
  
  
  
  #####get heatmap for protein units
  
  
  ptm_cox = fread ("ptm_cox_table.tsv")
  domain_cox = fread ("domain_cox_table.tsv")
  prot_cox = fread ("prot_cox_table.tsv")
  
  ptm_sig_cox <- ptm_cox %>% select(one_of(ptm_sig$identifier))
  
  trans_ptm = t(ptm_sig_cox[,-1])
  
  colnames(trans_ptm)= ptm_cox$patient_id
  
  trans_ptm_order = trans_ptm[, sel_patient]
  
  write.table(trans_ptm_order, "significant_ptm.tsv",quote= F, row.names = T,sep = "\t")
  
  
  #ptm_cb = c(0,0.5,1.5,2.5,3.5,4.5)
  
  
  
  #heatmap.2(trans_ptm_order, trace = "n", Rowv = T, Colv = F,col=colorRampPalette(c("white","black","yellow","orange","red")),breaks=ptm_cb,main = "ptm", dendrogram = "row",margins = c(2,14))
  
  
  ########
  
  
  #setwd("/Users/ginny/tcga_june_gastric_cancer/TCGA_STAD")
  domain_cox = fread ("domain_cox_table.tsv")
  prot_cox = fread ("prot_cox_table.tsv")
  
  domain_sig_cox <- domain_cox %>% select(one_of(domain_sig$identifier))
  
  trans_domain = t(domain_sig_cox[,-1])
  
  colnames(trans_domain)= domain_cox$patient_id
  
  trans_domain_order = trans_domain[, sel_patient]
  
  write.table(trans_domain_order, "significant_domain.tsv",quote= F, row.names =T,sep = "\t")
  
  
  mynames = rownames(trans_domain_order)
  get_prot = sapply(1:length(mynames),function(x) strsplit(mynames,split = "_")[[x]][1])
  
  write.table(unique(get_prot), "unique_significant_domain_names.tsv",quote= F, row.names = F,sep = "\t")
  
  write.table(unique(mynames), "unique_whole_significant_domain_names.tsv",quote= F, row.names = F,sep = "\t")
  
  
}




#domain_cb = c(0,0.5,1.5,2.5,3.5,4.5,5.5,6.5)
#heatmap.2(trans_domain_order, trace = "n", Rowv = T, Colv = F,col=colorRampPalette(c("white","black","yellow","orange","orange","red","red")),breaks=domain_cb,main = "domain", dendrogram = "row",margins = c(2,14))




#par(mfrow=(c(2,1)))

#heatmap.2(trans_ptm_order, trace = "n", Rowv = T, Colv = F,col=colorRampPalette(c("white","black","yellow","orange","red")),breaks=ptm_cb,main = "ptm", dendrogram = "row")

#heatmap.2(trans_domain_order, trace = "n", Rowv = T, Colv = F,col=colorRampPalette(c("white","black","yellow","orange","orange","red","red")),breaks=domain_cb,main = "domain", dendrogram = "row")


