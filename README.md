# Prepare tcga data

```{r}
source("prexa_tcgaPrep.R")

## need to direct to correct directory 
cancer_muse=fread("STAD_muse.maf", stringsAsFactors = F)
cancer_mutect=fread("STAD_mutect.maf", stringsAsFactors = F)
cancer_somatic=fread("STAD_somatic.maf", stringsAsFactors = F)
cancer_varscan=fread("STAD_varscan.maf", stringsAsFactors = F)
 
stomach_cancer = get_tcga_whole(cancer_muse, cancer_mutect, cancer_somatic, cancer_varscan, "STAD")
# 

```



# Map mutations to PTM sites
```{r}
source("prexa_ptmMapping.R")
ptm_data=fread("PSP_ptm_aggregate.tsv", header = T, stringsAsFactors = F)
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors= F)
get_ptm_mapping(stomach_all, "snp_ptm_mapping_sig_stad.tsv")
```


# Map mutations to domain
```{r}
source("prexa_domainMapping.R")
domain_data = fread ("domain_sorted_df.tsv", stringsAsFactors=F)
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors= F)
get_domain_mapping (stomach_all, "snp_domain_mapping_sig_stad.tsv")
```



# Parse patient info

```{r}

source("prexa_clinicalParsing.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
tcga_follow_up = fread("tcga_follow_up.txt", stringsAsFactors = F)
tcga_patient = fread("tcga_patient.txt", stringsAsFactors = F)
tcga_stad_patient_info = get_clinical(tcga_follow_up, tcga_patient, "stad_clinical.tsv")

```

# In case you have more than one follow up data, need to merge after dealing with each one

```{r}
source("prexa_clinicalMerging.R")
setwd("/data/ginny/PREXA/TCGA_BRCA")

mymerged = get_clinical_merged(c("brca_clinical_1_5.tsv",
                                  "brca_clinical_2_1.tsv", 
                                  "brca_clinical_4_0.tsv"),
                                "brca_clinical.tsv")
 
```





# Map domain mutations to patients
```{r}


source("prexa_domainPatientMapping.R")
setwd("/data/ginny/PREXA/TCGA_STAD")

sel_clinical = fread("stad_clinical.tsv", stringsAsFactors=F)

stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)

snp_domain_mapping = fread("snp_domain_mapping_sig.tsv", stringsAsFactors=F)
get_patient_domain_mapping(sel_clinical, 
stomach_all,
snp_domain_mapping,
"stad_patient_domain_table.tsv",
"stad_patient_domain_table_transpose.tsv")
```



# Map ptm mutations to patients
```{r}
source("prexa_ptmPatientMapping.R")
setwd("/data/ginny/PREXA/TCGA_STAD")

sel_clinical=fread("stad_clinical.tsv",stringsAsFactors = F)
stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)
snp_ptm_mapping=fread("snp_ptm_mapping_sig_stad.tsv", stringsAsFactors = F)

get_patient_ptm_mapping(sel_clinical, 
stomach_all,
snp_ptm_mapping,
"stad_patient_ptm_table.tsv",
"stad_patient_ptm_table_transpose.tsv")

```

# Map protein mutation to patients
```{r}
source("prexa_protPatientMapping.R")

setwd("/data/ginny/PREXA/TCGA_STAD")

sel_clinical = fread("stad_clinical.tsv", stringsAsFactors=F)

stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)
get_patient_prot_mapping(sel_clinical,
stomach_all,
"stad_patient_prot_table.tsv",
"stad_patient_prot_table_transpose.tsv")

```



# Build the patient_allLoci matrix

```{r}

source("prexa_wholeLoci.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors = F)
sel_clinical = fread ("stad_clinical.tsv",stringsAsFactors = F)

### focus only on mutations which happen on 10 or more patients

sel_loci <- stomach_all %>% group_by(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos,mutationaa, aapos, protlen)%>%filter(n()>=10)

###fix some possible bugs
#this_locus_patient <- sel_loci %>% mutate(locus_name = paste(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos, mutationaa,aapos,protlen,sep="-")) %>% ungroup() %>% select(locus_name, sampleid)


whole_stad_matrix = get_locus_patient_matrix(sel_loci,"stad",sel_clinical)



```

# Get selected cox matrix (exclude hypermutated patients)

```{r}
source("prexa_selectForCox.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
stad_whole = fread("stad_full_cox_table.tsv")
stad_clinical = fread("stad_clinical.tsv",stringsAsFactors = F)
get_selected_cox_matrix(stad_whole,  stad_clinical, 3500, 1, "stad")

```



#  SNP cox functions
```{r}

source("prexa_snpCox.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors = F)
sel_clinical = fread("stad_sel_sort_clinical.tsv", stringsAsFactors = F)

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


silent_coef_p_q = get_q_value("Silent", sel_loci,sel_clinical)
misSense_coef_p_q = get_q_value("Missense_Mutation",sel_loci,sel_clinical)
nonSense_coef_p_q = get_q_value("Nonsense_Mutation", sel_loci,sel_clinical)
spliceSite_coef_p_q= get_q_value("Splice_Site", sel_loci,sel_clinical)
inFrameDel_coef_p_q = get_q_value("In_Frame_Del", sel_loci,sel_clinical)
inFrameIns_coef_p_q = get_q_value("In_Frame_Ins", sel_loci,sel_clinical)
frameShiftDel_coef_p_q = get_q_value("Frame_Shift_Del", sel_loci,sel_clinical)
frameShiftIns_coef_p_q = get_q_value("Frame_Shift_Ins", sel_loci,sel_clinical)
nonStop_coef_p_q = get_q_value("Nonstop_Mutation", sel_loci,sel_clinical)
tFlank_coef_p_q = get_q_value("3'Flank", sel_loci,sel_clinical)
fFlank_coef_p_q = get_q_value("5'Flank", sel_loci,sel_clinical)
tUTR_coef_p_q = get_q_value("3'UTR", sel_loci,sel_clinical)
fUTR_coef_p_q = get_q_value("5'UTR", sel_loci,sel_clinical)
IGR_coef_p_q = get_q_value("IGR", sel_loci,sel_clinical)
intron_coef_p_q = get_q_value("Intron", sel_loci,sel_clinical)
RNA_coef_p_q = get_q_value("RNA", sel_loci,sel_clinical)
targetedRegion_coef_p_q = get_q_value("Targeted_Region", sel_loci,sel_clinical)
translationStartSite_coef_p_q = get_q_value("Translation_Start_Site", sel_loci,sel_clinical)


```


# Type(ptm, domain, prot) cox functions
```{r}

source("prexa_typeCox.R")
setwd("/data/ginny/PREXA/TCGA_STAD")


sel_clinical = fread("stad_sel_sort_clinical.tsv",stringsAsFactors = F)

sel_patients = sel_clinical$bcr_patient_barcode

df_sel_patients <- data.frame(sel_clinical.bcr_patient_barcode = sel_patients) %>% mutate_if (is.factor, as.character)





stad_patient_ptm = fread("stad_patient_ptm_table.tsv", stringsAsFactors = F, header = T)


sel_from_ptm = which(stad_patient_ptm$sel_clinical.bcr_patient_barcode %in% sel_patients)

stad_sel_patient_ptm = stad_patient_ptm [sel_from_ptm,]


sort_stad_sel_patient_ptm <- stad_sel_patient_ptm %>%
  left_join(df_sel_patients,.)
  
#######

stad_patient_domain =fread("stad_patient_domain_table.tsv", stringsAsFactors = F, header = T)

sel_from_domain = which(stad_patient_domain$sel_clinical.bcr_patient_barcode %in% sel_patients)

stad_sel_patient_domain = stad_patient_domain [sel_from_domain,]

sort_stad_sel_patient_domain <- stad_sel_patient_domain %>%
  left_join(df_sel_patients,.)
  

########
stad_patient_prot =fread("stad_patient_prot_table.tsv", stringsAsFactors = F, header = T)


sel_from_prot = which(stad_patient_prot$sel_clinical.bcr_patient_barcode %in% sel_patients)

stad_sel_patient_prot = stad_patient_prot [sel_from_prot,]


sort_stad_sel_patient_prot <- stad_sel_patient_prot %>%
  left_join(df_sel_patients,.)
  


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

time = sel_clinical$survival_time
status = sel_clinical$survival_ind


ptm_sel_coef_p_q = get_type_q_value(sort_stad_sel_patient_ptm,"ptm")
domain_sel_coef_p_q = get_type_q_value(sort_stad_sel_patient_domain,"domain")
prot_sel_coef_p_q = get_type_q_value(sort_stad_sel_patient_prot,"prot")




```

# Get the data for boxplot

```{r}

source("prexa_boxplotData.R")
setwd("/data/ginny/PREXA/TCGA_STAD")



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


stad_cox = get_df_all(silent_coef_p_q,
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




### this plot can be draw in my local R studio


stad_cox = fread("sel_all_coef_p_q.tsv", stringsAsFactors = F)


library(ggplot2)


ggplot(data = stad_cox, aes(x = all_labels, y = -log10(qValue)))+
  geom_boxplot(size = 0.3, outlier.size = 0.3)+
  ggtitle ("boxplot of -log10(qValue)")


```



# Get heatmap of signification genes

```{r}

source("prexa_sigHeatmapData.R")
setwd("/data/ginny/PREXA/TCGA_STAD")

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


get_data_for_heatmap (silent_coef_p_q,
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

```





