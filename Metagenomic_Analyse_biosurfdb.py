#!/usr/bin/env python
# coding: utf-8



# Pré-processamento dos dados
    # Import de biblio e leitura dos dados
import pandas as pd
dfs = pd.read_excel("/home/eulle/Documents/Master_Degree/resultados_biosurfdb_isolados_consorcios.xlsx",sheet_name=None)



# Criando uma lista de dataframes
df_sheets = []

for ii in dfs.keys():
    df_sheets.append(dfs[ii])

#def processing_data(dataframe):
df_sheets = [jj.rename(columns={"query (e.g., gene) sequence id": "query_sequence_id", ' \t subject (e.g., reference genome) sequence id': "subject_sequence_id", ' \t percentage of identical matches':'percentage of identical matches'})for jj in df_sheets]    
        
    

## PRÉ-PROCESSANDO OS DADOS(DADOS DE INTERESSE)
# NODE_1174_length_214_cov_1.65409 - bacillomycin_D_synthetase_A ( SÓ MUDA O QUERY START)
#Linhas com melhor e-value e/ou identidade
labels = df_sheets[1].query_sequence_id.unique()
df1 = pd.DataFrame(columns=df_sheets[1].columns)
df4 = pd.DataFrame(columns=df_sheets[1].columns)
for label in labels:
    subject_genes = df_sheets[1].loc[(df_sheets[1].query_sequence_id==label)]
    subject_genes = subject_genes['subject_sequence_id'].unique()
    subject_genes = [x.split("|")for x in subject_genes]
    subject_genes = list(set([ii[0]for ii in subject_genes]))
    for subject_gene in subject_genes:
        rows_subject_gene = df_sheets[1].loc[(df_sheets[1].query_sequence_id==label)]
        rows_subject_gene = rows_subject_gene[rows_subject_gene['subject_sequence_id'].str.contains(subject_gene+'+\|',regex=True)]

        value_df = rows_subject_gene['expect value'].min()
        if value_df == 0:
            df_evalue_0 = rows_subject_gene[rows_subject_gene['expect value']==0]
            df2 = df_evalue_0.loc[df_evalue_0.loc[:,'percentage of identical matches'].idxmax(),]
            df1 = df1.append(df2,ignore_index = True)

        else:
            df_evalue_diff = rows_subject_gene.loc[(rows_subject_gene.query_sequence_id==label)]
            df3 = df_evalue_diff.loc[df_evalue_diff.loc[:,'expect value'].idxmin(),]
            df4 = df4.append(df3,ignore_index=True)

df1 = df1.append(df4,ignore_index=True)
df1.to_csv('isolados_bg1.csv',index=False)


# GERAR AS IMAGENS
 

from dna_features_viewer import GraphicFeature, GraphicRecord
import random
import numpy as np

features = []
nodes = df1.query_sequence_id.unique()
#ff = df1[(df1.query_sequence_id=='NODE_116_length_960_cov_1.82541')&(df1['end of alignment in query'])].max()
for node in nodes:
    rows_nodes = df1[(df1.query_sequence_id==node)]
    for row in rows_nodes.iterrows():
        row[1]['subject_sequence_id'] = row[1]['subject_sequence_id'].split('|').pop(0)
    #[print(ii)for ii in row[1]['subject_sequence_id']]   
        features.append(GraphicFeature(start=row[1]['start of alignment in query'], end=row[1]['end of alignment in query'] , strand= +1 if row[1]['end of alignment in query'] > row[1]['start of alignment in query'] else -1 , color=list(np.random.choice(np.arange(0.0, 1.0, 0.1), size=3)),label=row[1]['subject_sequence_id']))
    record = GraphicRecord(sequence_length=rows_nodes['start of alignment in query'].max() + 100 if rows_nodes['start of alignment in query'].max() > rows_nodes['end of alignment in query'].max() else rows_nodes['end of alignment in query'].max() + 100, features=features)
    ax,_ = record.plot(figure_width=15)
    ax.figure.savefig(node+'.png')
    features = []


