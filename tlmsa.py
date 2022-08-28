'''
@author: onurdogan
'''

from bioservices import *
import pandas as pd
import numpy as np
from sumonet.utils.encodings import Encoding
from sumonet.model.architecture import SUMOnet


def getseq(df,uniprotCol):
    df['seq']= np.nan

    u = UniProt()
    unique_uniprotid=df[uniprotCol].unique()
    print(len(unique_uniprotid),'unique uniprotID')
    print('Started retrieving sequence from Uniprot')
    counter=0
    for j in unique_uniprotid:
        df_nn=df.loc[df[uniprotCol]==j]
        sequence2 = u.search(str(j), columns="sequence")
        for i in range(len(df_nn)):
            df.loc[df[uniprotCol]==j,('seq')]=str(sequence2[9:-1])
        counter=counter+1
        if (counter%100==0):
            print(len(unique_uniprotid)-counter,'left')
    return df

def getMutatedseq(df,caseIdCol,HugoSymbolCol,PositionCol,seqCol,aaNameCol):
    df['new_seq']=np.nan
    unique_case_id=df[caseIdCol].unique()
    for i in unique_case_id:
        df_n=df.loc[df[caseIdCol]==i]
        unique_genes_case_id=df_n[HugoSymbolCol].unique()
        for j in unique_genes_case_id:
            df_nn=df_n.loc[df[HugoSymbolCol]==j]
            df_nn=df_nn.reset_index(drop=True)
            
            seq=list(df_nn[seqCol][0])
            
            for k in range(int(len(df_nn))):
                n_pos=int(int(df_nn[PositionCol][k]))
                seq[int(n_pos)-1]=df_nn[aaNameCol][k]
            
            new_seq="".join(seq)
            
            for s in range(int(len(df_nn))):
                df.loc[(df[caseIdCol]==i) & (df[HugoSymbolCol]==j),'new_seq']=new_seq
    
    return df

def getSubseq(df,aaNameCol,PositionCol,newSeqCol):
    df['subseq']=np.nan
    count_row=df.shape[0]
    for i in range(count_row):
        df_n=df.iloc[[i]]
        if df_n[aaNameCol][i]=='K':
            
            pos=df_n[PositionCol]
            seq=str(df_n[newSeqCol].values)[2:-2]
            len_seq=len(seq)
            if int(pos)<11:
                a=seq[0:int(pos)+10]
                for s in range(int(11-pos)):
                    a="X"+a
                    
            elif int((len_seq-int(pos)))<11:
                a=seq[int(pos)-11:]
                
                for s in range(int(10-(len_seq-int(pos)))):
                    a=a+"X"
                
            else:
                a=seq[int(pos)-11:int(pos)+10]
            df.loc[int(i),'subseq']=str(a)
        else:
            continue
    return df


def motif_discover(motif):
    
    motif_inds = []
    motif_check = False
    y1 = ['I','L','V']
    y2 = ['A','F','I','L','M','P','V','W']
    y3 = ['A','F','G','I','L','M','P','V','W','Y']
    y4 = ['A','F','G','I','L','P','V']
    alpha = ['D','E']
    pS = ['S','T']
    
    if motif[9] in y1 and motif[12] in alpha: #motif1
        motif_inds.append(1) 
        motif_check = True
        
    if motif[9] in y2 and motif[12] in alpha:#motif2
        motif_inds.append(2) 
        motif_check = True
        
    if motif[9] in y3 and motif[12] in alpha:#motif3
        motif_inds.append(3) 
        motif_check = True
        
        
        
    if motif[9] in y2 and motif[12] in alpha and motif[15] == 'S' and motif[16] == 'P' :#motif4
        motif_inds.append(4) 
        motif_check = True
        
        
        
    if motif[9] in y2 and motif[12] in alpha: #motif5
        specific_motif = motif[14:-1]
        motif_counter = 0
        for aa in specific_motif:
            if aa in alpha:
                motif_counter += 1
        if motif_counter >= 2:
            motif_inds.append(5) 
            motif_check = True
        
        
        
    if motif[7] in y4 and motif[8] in y4 and motif[9] in y4 and motif[12] == 'E':#motif6
        motif_inds.append(6) 
        motif_check = True
        
        
        
    if motif[8] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[13] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[8] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[14] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[8] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[15] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[8] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[16] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[7] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[13] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[7] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[14] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[7] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[15] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[7] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[16] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[6] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[13] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[6] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[14] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[6] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[15] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[6] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[16] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[5] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[13] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[5] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[14] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[5] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[15] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[5] in ['P','G'] and motif[9] in ['I','V'] and motif[12]=='E' and motif[16] in ['P','G']:#motif7
        motif_inds.append(7) 
        motif_check = True
        
    if motif[9] in ['I','V'] and motif[12]=='E' and motif[13] == 'P':#motif8 -> expand
        motif_inds.append(8) 
        motif_check = True
        
    if motif[9] in ['I','V'] and motif[12]=='E' and motif[14] == 'P':#motif8 -> expand
        motif_inds.append(8) 
        motif_check = True
        
    if motif[9] in ['I','V'] and motif[12]=='E' and motif[15] == 'P':#motif8 -> expand
        motif_inds.append(8) 
        motif_check = True
        
    if motif[9] in ['I','V'] and motif[12]=='E' and motif[16] == 'P':#motif8 -> expand
        motif_inds.append(8) 
        motif_check = True
        
    if motif[9] in y2 and motif[12] in alpha and motif[13] == 'P':#motif9
        motif_inds.append(9) 
        motif_check = True
    
    if motif[9] in y2 and motif[12] == 'S' and motif[13] == 'P':#motif10
        motif_inds.append(10) 
        motif_check = True
        
    if motif[8] in alpha and motif[11] in y1:#motif11
        motif_inds.append(11) 
        motif_check = True
        
    if motif[8] in alpha and motif[11] in y2:#motif12
        motif_inds.append(12) 
        motif_check = True
        
    if motif[8] in alpha and motif[11] in y3:#motif13
        motif_inds.append(13) 
        motif_check = True
        
    if motif_check == False:
        motif_inds.append(0)
        
    return motif_inds

def motif_predict(seq):
    
    motif_names = ['non_motif','motif1','motif2','motif3','motif4','motif5','motif6',
               'motif7','motif8','motif9','motif10','motif11','motif12','motif13']    
    
    motifIndices = motif_discover(seq)
    
    if 0 not in motifIndices: #non motif 
    
        motifs = list(motif_names[i] for i in motifIndices)

        return (seq,motifs,1)
    
    else:
        
        return (seq,['non_motif'],0)

def get_motif(df,subseqCol):
    df['motif']=np.nan
    df['motif_names']=np.nan
    for i in range(df.shape[0]):
        a,b,c= motif_predict(df[subseqCol][i])
        df['motif'][i]=c
        df['motif_names'][i]=b
    return df
