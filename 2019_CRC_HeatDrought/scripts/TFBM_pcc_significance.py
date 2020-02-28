"""
Given a top hits file, add a column that says what level the PCC is significant to
based on the 95th percentile of PCCs between fam, within fam, and between random kmers.
"""
import pandas as pd
import numpy as np
import sys, os
from collections import defaultdict


FAM_KEY = '/mnt/home/azodichr/GitHub/MotifDiscovery/TF_Family_dictionary.txt'
EXP = '/mnt/home/azodichr/GitHub/MotifDiscovery/TFBM_with_bet_rand.txt'

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-top":
    TOP = sys.argv[i+1]  
  if sys.argv[i] == "-key":
    FAM_KEY = sys.argv[i+1]
  if sys.argv[i] == "-exp":    # Motifs to search through (subset of the -t2 tamo file)
    EXP = sys.argv[i+1]    
print(TOP)

if len(sys.argv) <= 1:
  print(__doc__)
  exit()

# Key to group TF families 
families = {}
with open(FAM_KEY,'r') as fam_file:
  for fline in fam_file:
    name = fline.strip().split('\t')[0]
    families[name] = fline.strip().split('\t')[1:]

# 
percen = pd.read_csv(EXP, header=0, sep='\t',index_col=None)
percen = pd.pivot_table(percen, values='95perc',index=['family','comparison'])
#print(test)

top = pd.read_csv(TOP, header=0, sep='\t', index_col=0)

sig = []
std_TF_name = []
for row in top.index:
  try: 
    PCC = top.at[row,'top_PCC']
    hit_fam = top.at[row,'top_hit']
  except:
    try:
      PCC = top.at[row,'DAP_top_PCC']
      hit_fam = top.at[row,'DAP_top_Fam']
    except:
      PCC = top.at[row,'CIS_top_PCC']
      hit_fam = top.at[row,'CIS_top_Fam']
  
  for TF_fam in families:
    if hit_fam in families[TF_fam]:
      hit_fam_standardized = TF_fam

  std_TF_name.append(hit_fam_standardized)
  if PCC >= percen[hit_fam_standardized]['within']:
    sig.append('within')
  elif PCC >= percen[hit_fam_standardized]['between']:
    sig.append('between')
  elif PCC >= percen[hit_fam_standardized]['random']:
    sig.append('random')
  else:
    sig.append('worse')
  
top['significance'] = sig
top['TF_Fam_standardized'] = std_TF_name
top.to_csv(TOP+"_sig", sep='\t')


print(top['significance'].value_counts())


#print(top.head())

