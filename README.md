# PREXA_Rversion
R code for my second project/part_1



## The pipeline of executing the scripts 


### 1 get the PTM and domain data of the proteins from previous analysis 
* prexa_ptmCollection.R
* prexa_domainCollection.R
### 2 map the TCGA SNP variance data 
* prexa_tcgaCollection.R
### 3 map TCGA variance to PTM and to domain on the proteins
* prexa_ptmMapping.R
* prexa_domainMapping.R


```{r}
source("prexa_ptmMapping.R")
ptm_data=fread("PSP_ptm_aggregate.tsv", header = T, stringsAsFactors = F)
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors= F)
get_ptm_mapping(stomach_all, "snp_ptm_mapping_sig_stad.tsv")
```
```{r}
source("prexa_domainMapping.R")
domain_data = fread ("domain_sorted_df.tsv", stringsAsFactors=F)
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors= F)
get_domain_mapping (stomach_all, "snp_domain_mapping_sig.tsv")
```



### 3 parse the clinical data from TCGA 
* prexa_clinicalParsing.R
```{r}

source("prexa_clinicalParsing.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
tcga_follow_up = fread("tcga_follow_up.txt", stringsAsFactors = F)
tcga_patient = fread("tcga_patient.txt", stringsAsFactors = F)
tcga_stad_patient_info = get_clinical(tcga_follow_up, tcga_patient, "stad_clinical.tsv")

```



### 4 map the SNP variance to both PTM/domain/protein and patient 
* prexa_ptmPatientMapping.R
* prexa_domainPatientMapping.R
* prexa_protPatientMapping.R



```{r}
source("prexa_ptmPatientMapping.R")
setwd("/data/ginny/PREXA/TCGA_STAD")

sel_clinical=fread("stad_clinical.tsv",stringsAsFactors = F)
stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)
snp_ptm_mapping=fread("snp_ptm_mapping_sig_stad.tsv", stringsAsFactors = F)

get_patient_ptm_mapping(sel_clinical, stomach_all, snp_ptm_mapping,
"stad_patient_ptm_table.tsv",
"stad_patient_ptm_table_transpose.tsv")

```
```{r}


source("prexa_domainPatientMapping.R")
setwd("/data/ginny/PREXA/TCGA_STAD")

sel_clinical = fread("stad_clinical.tsv", stringsAsFactors=F)

stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)

snp_domain_mapping = fread("snp_domain_mapping_sig.tsv", stringsAsFactors=F)
get_patient_domain_mapping(sel_clinical, stomach_all, snp_domain_mapping,
"stad_patient_domain_table.tsv","stad_patient_domain_table_transpose.tsv")
```

```{r}
source("prexa_protPatientMapping.R")

setwd("/data/ginny/PREXA/TCGA_STAD")

sel_clinical = fread("stad_clinical.tsv", stringsAsFactors=F)

stomach_all=fread("STAD_whole.tsv",stringsAsFactors = F)
get_patient_prot_mapping(sel_clinical, stomach_all,
"stad_patient_prot_table.tsv","stad_patient_prot_table_transpose.tsv")

```

### 5 build a matrix where each row is a mutation and each column is a patient

```{r}

source("prexa_wholeLoci.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
stomach_all = fread("STAD_whole.tsv", stringsAsFactors = F)
sel_clinical = fread ("stad_clinical.tsv",stringsAsFactors = F)

### focus only on mutations which happen on 10 or more patients

sel_loci <- stomach_all %>% group_by(genename, uni_swiss, uni_trembl,chromosome, mutationtype, mutationpos,mutationaa, aapos, protlen)%>%filter(n()>=10)

whole_stad_matrix = get_locus_patient_matrix(sel_loci,"stad",sel_clinical)

```

### 6 get selected cox matrix and sorted patients

```{r}
source("prexa_selectForCox.R")
setwd("/data/ginny/PREXA/TCGA_STAD")
stad_whole = fread("stad_full_cox_table.tsv")
stad_clinical = fread("stad_clinical.tsv",stringsAsFactors = F)
get_selected_cox_matrix(stad_whole,  stad_clinical, 3500, 1, "stad")

```


### 7 perform cox regression with the data matrix acquired on pure SNP data
* prexa_snpCox.R
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
```



### 8 perform cox regression with the data matrix acquired from PTM/domain/protein data
* prexa_typeCox.R

```{r}

source("prexa_typeCox.R")
setwd("/data/ginny/PREXA/TCGA_STAD")

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

time = sel_clinical$survival_time
status = sel_clinical$survival_ind



ptm_coef_p_q = get_type_q_value(stad_patient_ptm,"ptm")
domain_coef_p_q = get_type_q_value(stad_patient_domain,"domain")
prot_coef_p_q = get_type_q_value(stad_patient_prot,"prot")


```


### 9 visualize with boxplot
* prexa_boxplot.R
