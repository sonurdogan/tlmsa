from bioservices import *
import pandas as pd
import numpy as np
import tlmsa
from sumonet.utils.encodings import Encoding
from sumonet.model.architecture import SUMOnet


df=pd.read_csv('k_all_mutation2.csv')
df=df.drop_duplicates(subset=['case_id', 'Hugo_Symbol','HGVSp_Short'])
df=df.reset_index(drop=True)

df=df[:300]


df2=tlmsa.getseq(df,'uniprotID')
df3=tlmsa.getMutatedseq(df2,'case_id','Hugo_Symbol','positions_','seq','aa_name')
df4=tlmsa.getSubseq(df3,'aa_name','positions_','new_seq')

df4=df4[['case_id','Hugo_Symbol','Tumor_Sample_Barcode','positions_','aa_name','subseq']]
df4.to_csv('mutatedseq_res.csv',index=False)


df_sumo=df4[df4.aa_name=='K']
df_sumo=df_sumo.reset_index(drop=True)
df_sumo=tlmsa.predict_sumo(df_sumo,'subseq')
df_res=tlmsa.predict_motif(df_sumo,'subseq')
df_res=df_res.sort_values(by=['Sumo_prob'],ascending=False)

df_res.to_csv('sumo_result.csv',index=False)

