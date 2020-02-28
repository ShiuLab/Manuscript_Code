""" 
PURPOSE:
Find best match TFBM from PCC-Distance Matrix or matricies


INPUT:
  -df1      
  -df2        
  -save  

OUTPUT:
  -save.txt

AUTHOR: Christina Azodi

REVISIONS:   Submitted 8/24/2017
"""

import pandas as pd
import numpy as np
import sys, os

CIS_KEY = '/mnt/home/mjliu/kmer_5/Athaliana_TFBM_v1.01.tm.index.direct.index'
DAP_KEY = '/mnt/research/ShiuLab/14_DAPseq/PWM_to_tamo/DAP_motifs.txt.tm_index'
KEEP = 'all'
DAP = CIS = 'skip'

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-t":
    T = sys.argv[i+1]  
  if sys.argv[i] == "-pcc_dap":
    DAP = sys.argv[i+1]
  if sys.argv[i] == "-pcc_cis":
    CIS = sys.argv[i+1]
  if sys.argv[i] == "-cis_key":   # Optional - default above
    CIS_KEY = sys.argv[i+1]
  if sys.argv[i] == "-keep":    # Motifs to search through (subset of the -t2 tamo file)
    KEEP = sys.argv[i+1]    
  if sys.argv[i] == '-save':
    SAVE = sys.argv[i+1]

if len(sys.argv) <= 1:
  print(__doc__)
  exit()

# Get column names for DAP and CIS-BP
dap_names = open(DAP_KEY,'r').readlines()
dap_names = [item.strip().split('\t')[0] for item in dap_names]
dap_names.remove('#motif_id')
cis_names = open(CIS_KEY,'r').readlines()
cis_names = [item.strip().split('\t')[0] for item in cis_names]

# Get row names for whatever you mapped
index = []
with open(T, 'r') as tamo_file:
  for l in tamo_file.readlines():
    if l.startswith('Log-odds matrix for Motif'):
      index.append(l.strip().split(' ')[-2])


if KEEP != 'all':
  keep = open(KEEP, 'r').readlines()
  keep = [item.strip() for item in keep]


if DAP != 'skip':
  dap = pd.read_csv(DAP, header=None, sep='\t', index_col=False, names=dap_names)
  dap.index = index
  #save_nam_dap = str(DAP) + '_labeled'
  #dap.to_csv(save_nam_dap, sep='\t')
  
  if KEEP != 'all':
    dap_keep_list = []
    for col in dap.columns:
      name = col.split('_')[1]
      if name in keep:
        dap_keep_list.append(col)
    dap = dap[dap_keep_list]

  # Get top DAP hit and calculate PCC from PCCD
  df = pd.DataFrame(0, index=dap.index.tolist(), columns = ['DAP_top'])

  if len(dap.columns) ==0:
    DAP='skip'

  elif len(dap.columns) == 1:
    df['DAP_top'] = dap.columns[0]
    df['DAP_top_PCCD'] = dap[dap.columns[0]]
    df['DAP_top_PCC'] = 1 - df['DAP_top_PCCD']
     
  elif len(dap.columns) > 1:
    df['DAP_top'] = dap.idxmin(axis=1)
    df['DAP_top_PCCD'] = dap.min(axis=1)
    df['DAP_top_PCC'] = 1 - df['DAP_top_PCCD'] 
    df['DAP_top_Fam'], df['DAP_TF'] = df['DAP_top'].str.split('_',1).str



if CIS != 'skip':
  cis = pd.read_csv(CIS, header=None, sep='\t', index_col=False, names = cis_names)
  print(cis.head())
  cis.index = index
  #save_nam_cis = str(CIS) + '_labeled'
  #cis.to_csv(save_nam_cis, sep='\t')
  
  if KEEP != 'all':
    cis_keep_list = []
    for col in cis.columns:
      name = col.split('.')[0] + '.02'
      if name in keep:
        cis_keep_list.append(col)
    cis = cis[cis_keep_list]
  cis_key = pd.read_csv(CIS_KEY, header=None, sep='\t', index_col=None)
  cis_key.columns = ['ID','PWM','Genes','CIS_Fam']

  if DAP == 'skip': #If df isn't already created by DAP:
    df = pd.DataFrame(0, index=cis.index.tolist(), columns = ['CIS_top'])

  if len(cis.columns) == 0:
    CIS='skip'

  elif len(cis.columns) == 1:
    df['CISBP_top'] = cis.columns[0]
    df['CISBP_top_PCCD'] = cis[cis.columns[0]]
    df['CISBP_top_PCC'] = 1 - df['CISBP_top_PCCD']
    df = df.merge(cis_key, how = 'left', left_on = 'CISBP_top', right_on = 'ID', left_index=True)

  elif len(cis.columns) > 1:
    df['CISBP_top'] = cis.idxmin(axis=1)
    df['CISBP_top_PCCD'] = cis.min(axis=1)
    df['CISBP_top_PCC'] = 1 - df['CISBP_top_PCCD']
    df = df.merge(cis_key, how = 'left', left_on = 'CISBP_top', right_on = 'ID', left_index=True)
  

# Reset the df index to the pCREs
try: 
  df['pCRE'] = dap.index.tolist()
except: 
  df['pCRE'] = cis.index.tolist()
df = df.set_index(keys ='pCRE', drop=True)

if DAP != 'skip' and CIS != 'skip':
  df['top_hit'] = np.where((df['DAP_top_PCC'] >= df['CISBP_top_PCC']),
    df['DAP_top_Fam'],df['CIS_Fam'])
  df['top_PCC'] = np.where((df['DAP_top_PCC'] >= df['CISBP_top_PCC']),
    df['DAP_top_PCC'],df['CISBP_top_PCC'])
  df['top_type'] = np.where((df['DAP_top_PCC'] >= df['CISBP_top_PCC']),
    'DAP','CISBP')

print(df.head())
save_name = SAVE + '_TopHits.txt'
df.to_csv(save_name,sep='\t')


