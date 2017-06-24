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
### 3 parse the clinical data from TCGA 
* prexa_clinicalParsing.R
### 4 map the SNP variance to both PTM/domain/protein and patient 
* prexa_ptmPatientMapping.R
* prexa_domainPatientMapping.R
### 5 perform cox regression with the data matrix acquired on pure SNP data
* prexa_snpCox.R
### 6 perform cox regression with the data matrix acquired from PTM/domain/protein data
* prexa_typeCox.R
### 7 visualize with boxplot
* prexa_boxplot.R
