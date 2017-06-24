## 3 processing the TCGA_SNP data


library(dplyr)
library(magrittr)
library(data.table)



set_whole_data=function(tcga)
{
    genename=tcga$Hugo_Symbol
    uni_swiss=tcga$SWISSPROT
    uni_trembl=tcga$TREMBL
    chromosome = tcga$Chromosome
    mutationtype=tcga$Variant_Classification
    mutationpos=tcga$Start_Position
    protpos=tcga$Protein_position
    mutationaa=tcga$Amino_acids
    samplecode=tcga$Tumor_Sample_Barcode
    
    getpatient=function(samplecodeone)
    {
        spt=unlist(strsplit(samplecodeone,split="-"))
        patientid=paste(spt[1],spt[2],spt[3],sep="-")
        patientid
    }
    
    
    get_tumor_status=function(samplecodeone)
    {
        spt=unlist(strsplit(samplecodeone,split="-"))
        status=spt[4]
        status
    }
    
    sampleid=unlist(lapply(samplecode, getpatient))
    sample_status=unlist(lapply(samplecode,get_tumor_status))
    
    df=data.frame(genename,uni_swiss,uni_trembl,chromosome,mutationtype,mutationpos,protpos,mutationaa,sampleid,sample_status, stringsAsFactors = F)
    
    bothin=strsplit(as.character(df$protpos),split="/")
    aapos=unlist(lapply(bothin, function(x){x[1]}))
    protlen=unlist(lapply(bothin, function(x){x[2]}))
    
    newdf=data.frame(df,aapos,protlen,stringsAsFactors = F)
    
    newdf
    
}

get_tcga_whole = function(cancer_muse, cancer_mutect, cancer_somatic, cancer_varscan, cancer_type)
{
    
    cancer_muse_df=set_whole_data(cancer_muse)
    cancer_mutect_df=set_whole_data(cancer_mutect)
    cancer_somatic_df=set_whole_data(cancer_somatic)
    cancer_varscan_df=set_whole_data(cancer_varscan)
    
    cancer_snp_df=unique(rbind(cancer_muse_df,cancer_mutect_df, cancer_somatic_df, cancer_varscan_df))
    write.table(cancer_snp_df, paste0(cancer_type, "_whole.tsv"),quote= F, row.names = F, sep="\t")
    
    
}



#cancer_muse=fread("STAD_muse.maf", stringsAsFactors = F)
#cancer_mutect=fread("STAD_mutect.maf", stringsAsFactors = F)
#cancer_somatic=fread("STAD_somatic.maf", stringsAsFactors = F)
#cancer_varscan=fread("STAD_varscan.maf", stringsAsFactors = F)

#stomach_cancer = get_tcga_whole(cancer_muse, cancer_mutect, cancer_somatic, cancer_varscan, "STAD")
