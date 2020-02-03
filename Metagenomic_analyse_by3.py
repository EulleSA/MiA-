#!/usr/bin/env python
# coding: utf-8


# Pré-processamento dos dados
    # Import de biblio e leitura dos dados
import pandas as pd
dfs = pd.read_excel("/home/eulle/Documents/Master_Degree/resultados_dos_passos_bioinformatica_realizados_isolados_consorcio_by3.xlsx",sheet_name=None)

# O sheet de interesse é o 11 -- by3
df_sheets = []
for ii in dfs.keys():
    df_sheets.append(dfs[ii])
    
# alterar o nome da coluna
df_sheets = [jj.rename(columns={"#query_name": "query_name"})for jj in df_sheets]    


import re
# definir os genes sem dados para 'no_gene' e capturar os dados de inicio e fim dos genes    
df_sheets[11]['Preferred_name'] = df_sheets[11]['Preferred_name'].fillna('no_gene')
gene_pos = [re.findall(r"^.*_(.*)_(.*)_(.)$", ii[1]['query_name']) for ii in df_sheets[11].iterrows()]
df_sheets[11]['query_name'].replace(regex=True,inplace=True,to_replace=r"_\d+_\d+_(.)$",value=r'')
gene_pos = [list(gene_pos[jj][0]) for jj in range(0,len(gene_pos))]

# Reverter as posições que tem o sinal '-'
start = []
end = []
for pos in gene_pos:
    if pos[2] == '-':
        pos[0],pos[1] = pos[1], pos[0]
        
    start.append(pos[0])
    end.append(pos[1])
df_sheets[11]['start_query'] = list(map(int,start))
df_sheets[11]['end_query'] = list(map(int,end))



# GERAR O VISUALIZADOR DINÂMICO DE ORFs

from dna_features_viewer import GraphicFeature, GraphicRecord,CircularGraphicRecord
from bokeh.resources import CDN
from bokeh.embed import file_html

import random
import numpy as np


features = []
nodes = df_sheets[11].query_name.unique()


rows_nodes = df_sheets[11][(df_sheets[11].query_name=='NODE_1_length_601576_cov_48.1615')]
for row in rows_nodes.iterrows():
    features.append(GraphicFeature(start=row[1]['start_query'], end=row[1]['end_query'] , strand= +1 if row[1]['end_query'] > row[1]['start_query'] else -1 , color='#fffaaa',label=row[1]['Preferred_name']))
record = GraphicRecord(sequence_length=20000, features=features)
plot = record.plot_with_bokeh(figure_width=10)
with open("testa.html", "w+") as f:
    f.write(file_html(plot, CDN, "Example Sequence"))


