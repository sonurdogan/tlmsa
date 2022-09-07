import pandas as pd
import numpy as np
import tlmsa
from sumonet.utils.encodings import Encoding
from sumonet.model.architecture import SUMOnet


df=pd.read_csv('TCGA-STAD.csv')
df=df.drop_duplicates(subset=['case_id', 'Hugo_Symbol','HGVSp_Short'])
df=df.reset_index(drop=True)


df2=tlmsa.getseq(df,'uniprotID')
df3=tlmsa.getMutatedseq(df2,'case_id','Hugo_Symbol','positions_','seq','aa_name')
df4=tlmsa.getSubseq(df3,'aa_name','positions_','new_seq')
df4=df4.dropna(subset=['subseq'])
df4=df4.reset_index(drop=True)

df4=df4[['case_id','Hugo_Symbol','Tumor_Sample_Barcode','positions_','aa_name','subseq','seq','new_seq']]
df4.to_csv('STAD-seq.csv',index=False)



df4=df4[['case_id','Hugo_Symbol','Tumor_Sample_Barcode','positions_','aa_name','subseq']]
df_sumo=df4[df4.aa_name=='K']
df_sumo=df_sumo.reset_index(drop=True)
df=tlmsa.get_motif(df_sumo,'subseq')


#predict possible sumoylation using SUMOnet
x=df['subseq'].values.tolist()
s=len(x)
list_of_ones= [1] * s
y=list_of_ones

encoder = Encoding(encoderType='blosum62') 
x_Test,y_Test = encoder.get_encoded_vectors_from_data(x,y)
input_shape=x_Test.shape

SUMOnet3_model = SUMOnet()
SUMOnet3_model.build(input_shape)
SUMOnet3_model.load_weights()  
y_preds = SUMOnet3_model.predict(x_Test)

#add prediction probabilities to dataframe

df['nonSumo_prob']=np.nan
df['Sumo_prob']=np.nan 
for i in range(s):
    df['nonSumo_prob'][i]=y_preds[i][0]
    df['Sumo_prob'][i]=y_preds[i][1]


df_res=df.sort_values(by=['Sumo_prob'],ascending=False)

df_res.to_csv('STAD-result.csv',index=False)
