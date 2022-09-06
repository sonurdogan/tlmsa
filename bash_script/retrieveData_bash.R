library(TCGAbiolinks)
library(DT)
library(readr)
library(stringr)
library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

#define query
query <- GDCquery(
  project =args[2],
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation",
  sample.type="Primary Tumor")

#download data to local computer/ If downloaded, no need to download again.
GDCdownload(query)

#prepare/get GDC data
maf <- GDCprepare(query)

#get just following features
maf=maf[,c('case_id','Hugo_Symbol','Tumor_Sample_Barcode','Transcript_ID','Variant_Classification','HGVSp','HGVSp_Short','SWISSPROT')]

#get mutation to lysines
lysines=as.data.frame(maf[maf$HGVSp_Short %like% "K$",c('case_id','Hugo_Symbol','Tumor_Sample_Barcode','Transcript_ID','Variant_Classification','HGVSp','HGVSp_Short','SWISSPROT')])

#get all mutation of genes that has mutation to lysines
unique_case_id=unique(lysines$case_id)

mat = matrix(ncol = 8)
df=data.frame(mat)
colnames(df)=colnames(maf)


for (i in unique_case_id) {
  r= maf[maf$case_id==i,]
  s= unique(lysines[lysines$case_id==i,'Hugo_Symbol'])
  
  for (j in s){
    m=r[r$Hugo_Symbol==j,]
    
    df=rbind(df,m)
  }
}

#get just missense mutation ones
df=df[df$Variant_Classification=='Missense_Mutation',]
data=as.data.frame(df[,c('case_id','Hugo_Symbol','Tumor_Sample_Barcode','Transcript_ID','HGVSp_Short','SWISSPROT')])

#extract position of mutation
s_pos=substr(data$HGVSp_Short,3,30)
positions_=parse_number(s_pos)

#extract uniprotID
uniprotID=word(data$SWISSPROT,1,sep = "\\.")

#extract amino acid name
aa_name=str_sub(data$HGVSp_Short,-1,-1)

#combine new features
data_new=cbind(data,aa_name,positions_,uniprotID)

# drop the rows with NA value in uniprotID
data_new <- data_new %>% drop_na(uniprotID)

#get just following features all data
data_new=data_new[,c('case_id','Hugo_Symbol','Tumor_Sample_Barcode','Transcript_ID','HGVSp_Short','aa_name','positions_','uniprotID')]
#save the data
str2=".csv"
file_name = paste(args[2],str2,sep="")

write.csv(data_new, file = file_name,row.names=FALSE)
