library(TCGAbiolinks)
library(DT)
library(readr)
library(stringr)
library(tidyverse)
library(data.table)

#define query
query <- GDCquery(
  project = "TCGA-SARC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation",
  sample.type="Primary Tumor")

#download data to local computer/ If downloaded, no need to download again.
GDCdownload(query)

#prepare/get GDC data default
maf <- GDCprepare(query)

#get just following features
maf=maf[,c('case_id','Hugo_Symbol','SYMBOL','Gene','Tumor_Sample_Barcode','Transcript_ID','Variant_Classification','Variant_Type','Amino_acids','Protein_position','HGVSp_Short','SWISSPROT','UNIPROT_ISOFORM')]

#get mutation to lysines
lysines=as.data.frame(maf[maf$HGVSp_Short %like% "K$",c('case_id','Hugo_Symbol','SYMBOL','Gene','Tumor_Sample_Barcode','Transcript_ID','Variant_Classification','Variant_Type','Amino_acids','Protein_position','HGVSp_Short','SWISSPROT','UNIPROT_ISOFORM')])

#get all mutation of genes that has mutation to lysines
unique_case_id=unique(lysines$case_id)

mat = matrix(ncol = 13)
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

data=df
data=data[-1,]

#delete the cases that has these mutations
del_mutations=list("Splice_Region","Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")

new_mat=matrix(ncol=13)
df_new=data.frame(new_mat)
colnames(df_new)=colnames(data)
df_new=df_new[-1,]

for (i in unique(data$case_id)) {
  r= data[data$case_id==i,]
  s= unique(data[data$case_id==i,'Hugo_Symbol'])
  for (j in s){
    m=r[r$Hugo_Symbol==j,]
    mutations=unique(m$Variant_Classification)
    print(mutations)
    if(length(intersect(mutations, del_mutations)) == 0){
        df_new=rbind(df_new,m)
    }
    else{
      next
    }
  }
}

data=df_new
data<-data[(data$Variant_Classification=="In_Frame_Del" | data$Variant_Classification=="In_Frame_Ins" | data$Variant_Classification=="Missense_Mutation"),]


#extract position of mutation
s_pos=substr(data$HGVSp_Short,3,30)
Position=parse_number(s_pos)

#extract uniprotID
UniprotID=word(data$SWISSPROT,1,sep = "\\.")

#extract amino acid name
Emerged_aa=str_sub(data$HGVSp_Short,-1,-1)

#combine new features
data_new=cbind(data,Emerged_aa,Position,UniprotID)


# drop the rows with NA value in uniprotID
data_new <- data_new %>% drop_na(UniprotID)

#get just following features all data
data_new=data_new[,c('case_id','Hugo_Symbol','SYMBOL','Gene','Tumor_Sample_Barcode','Transcript_ID','Variant_Classification','Variant_Type','Amino_acids','Emerged_aa','Protein_position','Position','HGVSp_Short','UNIPROT_ISOFORM','UniprotID')]
#save the data
write.csv(data_new, file = "TCGA-SARC.csv",row.names=FALSE)

