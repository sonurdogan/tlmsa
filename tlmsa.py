'''
@author: onurdogan

'''

from bioservices import *
import pandas as pd
import numpy as np
from sumonet.utils.encodings import Encoding
from sumonet.model.architecture import SUMOnet

class TLMSA:
    def __init__(self, df):
        self.data = df
        self.new_seqCol = 'mutated_seq'
        self.subseqCol = 'subseq'
    def get_sequence(self, uniprotCol):
        self.data['seq'] = np.nan

        u = UniProt()
        unique_uniprotid = self.data[uniprotCol].unique()
        print(len(unique_uniprotid), 'unique uniprotID detected')
        print('Started retrieving sequence from Uniprot!')
        counter = 0
        for j in unique_uniprotid:
            df_nn = self.data.loc[self.data[uniprotCol] == j]
            sequence2 = u.search(str(j), columns="sequence")
            for i in range(len(df_nn)):
                self.data.loc[self.data[uniprotCol] == j, ('seq')] = str(sequence2[9:-1])
            counter = counter + 1
            if (counter % 100 == 0):
                print(len(unique_uniprotid) - counter, 'left')

    def get_mutated_sequence(self, caseIdCol, HugoSymbolCol, PositionCol, seqCol, aaNameCol,new_seqCol=''):
        self.data[new_seqCol] = np.nan
        if new_seqCol!='':
            self.new_seqCol=new_seqCol

        self.PositionCol=PositionCol
        self.aaNameCol=aaNameCol

        unique_case_id = self.data[caseIdCol].unique()
        for i in unique_case_id:
            df_n = self.data.loc[self.data[caseIdCol] == i]
            unique_genes_case_id = df_n[HugoSymbolCol].unique()
            for j in unique_genes_case_id:
                df_nn = df_n.loc[self.data[HugoSymbolCol] == j]
                df_nn = df_nn.reset_index(drop=True)

                seq = list(df_nn[seqCol][0])

                for k in range(int(len(df_nn))):
                    n_pos = int(int(df_nn[PositionCol][k]))
                    if n_pos > len(seq):
                        print('position is out of sequence')
                        continue
                    else:
                        seq[int(n_pos) - 1] = df_nn[aaNameCol][k]

                new_seq = "".join(seq)

                for s in range(int(len(df_nn))):
                    self.data.loc[(self.data[caseIdCol] == i) & (self.data[HugoSymbolCol] == j), self.new_seqCol] = new_seq

    def get_subsequence(self, aaNameCol='', PositionCol='',new_seqCol='',subseqCol=''):
        if aaNameCol!='':
             self.aaNameCol=aaNameCol
        if PositionCol!='':
            self.PositionCol=PositionCol
        if new_seqCol!='': 
            self.new_seqCol=new_seqCol
        if subseqCol!='': 
            self.subseqCol=subseqCol
        
        self.data[self.subseqCol] = np.nan
        count_row = self.data.shape[0]
        for i in range(count_row):
            df_n = self.data.iloc[[i]]
            if df_n[self.aaNameCol][i] == 'K':
                pos = df_n[self.PositionCol]
                seq = str(df_n[self.new_seqCol].values)[2:-2]
                len_seq = len(seq)
                if len_seq < int(pos):
                    continue
                else:
                    if int(pos) < 11:
                        a = seq[0:int(pos) + 10]
                        for s in range(int(11 - pos)):
                            a = "X" + a

                    elif int((len_seq - int(pos))) < 11:
                        a = seq[int(pos) - 11:]
                        for s in range(int(10 - (len_seq - int(pos)))):
                            a = a + "X"

                    else:
                        a = seq[int(pos) - 11:int(pos) + 10]
                    self.data.loc[int(i), self.subseqCol] = str(a)
            else:
                continue
    def remove_nonK_aa(self, aaNameCol=''):
        if aaNameCol!='': 
            self.aaNameCol=aaNameCol
        
        self.data=self.data[self.data[self.aaNameCol]=='K']
        self.data=self.data.reset_index(drop=True)

    def get_motif(self,subseqCol=''):
        if subseqCol!='': 
            self.subseqCol=subseqCol

        self.data['motif'] = np.nan
        self.data['motif_names'] = np.nan
        for i in range(self.data.shape[0]):
            a, b, c = self.motif_predict(self.data[self.subseqCol][i])
            self.data['motif'][i] = c
            self.data['motif_names'][i] = b
    
    def get_sumonet_res(self,subseqCol='',EncoderType='blosum62'):

        if subseqCol!='': 
            self.subseqCol=subseqCol

        x=self.data[self.subseqCol].values.tolist()

        encoder = Encoding(encoderType=EncoderType) 
        x_Test = encoder.encode_data(x)
        input_shape=x_Test.shape

        SUMOnet3_model = SUMOnet()
        SUMOnet3_model.build(input_shape)
        SUMOnet3_model.load_weights()  
        y_preds = SUMOnet3_model.predict(x_Test)
        self.data['nonSumo_prob']=np.nan
        self.data['Sumo_prob']=np.nan 

        for i in range(len(x)):
            self.data['nonSumo_prob'][i]=y_preds[i][0]
            self.data['Sumo_prob'][i]=y_preds[i][1]



    def motif_discover(self, motif):
        
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

    def motif_predict(self, seq):
    
        motif_names = ['non_motif','motif1','motif2','motif3','motif4','motif5','motif6',
               'motif7','motif8','motif9','motif10','motif11','motif12','motif13']    
    
        motifIndices = self.motif_discover(seq)
    
        if 0 not in motifIndices: #non motif 
    
            motifs = list(motif_names[i] for i in motifIndices)

            return (seq,motifs,1)

        else:
        
            return (seq,['non_motif'],0)