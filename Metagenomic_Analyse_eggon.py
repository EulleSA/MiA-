#!/usr/bin/env python
# coding: utf-8


# Pré-processamento dos dados
    # Import de biblio e leitura dos dados
import pandas as pd
dfs = pd.read_excel("/home/eulle/Documents/Master_Degree/resultados_dos_passos_bioinformatica_realizados_isolados_consorcio.xlsx",sheet_name=None)

# O sheet de interesse é o 8 -- bg1   
df_sheets = []
for ii in dfs.keys():
    df_sheets.append(dfs[ii])
    
# alterar o nome da coluna
df_sheets = [jj.rename(columns={"#query_name": "query_name"})for jj in df_sheets]    


import re    
    
# definir os genes sem dados para 'no_gene' e capturar os dados de inicio e fim dos genes    
df_sheets[8]['Preferred_name'] = df_sheets[8]['Preferred_name'].fillna('no_gene')
gene_pos = [re.findall(r"^.*_(.*)_(.*)_(.)$", ii[1]['query_name']) for ii in df_sheets[8].iterrows()]
df_sheets[8]['query_name'].replace(regex=True,inplace=True,to_replace=r"_\d+_\d+_(.)$",value=r'')
gene_pos = [list(gene_pos[jj][0]) for jj in range(0,len(gene_pos))]

# Reverter as posições que tem o sinal '-'
start = []
end = []
for pos in gene_pos:
    if pos[2] == '-':
        pos[0],pos[1] = pos[1], pos[0]
        
    start.append(pos[0])
    end.append(pos[1])
df_sheets[8]['start_query'] = list(map(int,start))
df_sheets[8]['end_query'] = list(map(int,end))


# Pegar os melhores genes
labels = df_sheets[8].query_name.unique()

df4 = pd.DataFrame(columns=df_sheets[8].columns)

for label in labels:
    subject_genes = df_sheets[8].loc[(df_sheets[8].query_name==label)]
    subject_genes = subject_genes['Preferred_name'].unique()
    for subject_gene in subject_genes:
        rows_subject_gene = df_sheets[8].loc[(df_sheets[8].query_name==label)]
        rows_subject_gene = rows_subject_gene[rows_subject_gene['Preferred_name'].str.contains(subject_gene,regex=True)]
        
        df_evalue_diff = rows_subject_gene.loc[(rows_subject_gene.query_name==label)]
        df3 = df_evalue_diff.loc[df_evalue_diff.loc[:,'seed_ortholog_evalue'].idxmin(),]
        df4 = df4.append(df3,ignore_index=True)
df4.to_csv('isolados_bg1_eggon.csv',index=False)


# GERAR AS IMAGENS

from dna_features_viewer import GraphicFeature, GraphicRecord
import random
import numpy as np


features = []
nodes = df4.query_name.unique()
for node in nodes:
    rows_nodes = df4[(df4.query_name==node)]
    for row in rows_nodes.iterrows():
        features.append(GraphicFeature(start=row[1]['start_query'], end=row[1]['end_query'] , strand= +1 if row[1]['end_query'] > row[1]['start_query'] else -1 , color=list(np.random.choice(np.arange(0.0, 1.0, 0.1), size=3)),label=row[1]['Preferred_name']))
    record = GraphicRecord(sequence_length=rows_nodes['start_query'].max() + 100 if rows_nodes['start_query'].max() > rows_nodes['end_query'].max() else rows_nodes['end_query'].max() + 100, features=features)
    ax,_ = record.plot(figure_width=15)
    ax.figure.savefig('plot_eggon/'+node+'.png')
    features = []

