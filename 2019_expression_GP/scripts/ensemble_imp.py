"""
Take multiple _imp and _coef.csv files and 
generate an ensemble importance score.

Approach:
  1. Calculate mean score for each model (already done for _imp) 
  2. Normalize those scores between 0 and 1
  3. Take the mean of the normalized scores

python ensemble_imp.py SAVE_NAME FILES...

After SAVE_NAME arg, can give as many _imp and _coef.csv files as desired

"""
import sys, os
import pandas as pd
import numpy as np

n = 100
SAVE = sys.argv[1]

# Read in scores from models to ensemble
items = []
for i in range (2,len(sys.argv)):
  items.append(sys.argv[i])

print(items)


# For each replicate

count = 1
for f in items:
  print(f)
  if "_coef.csv" in f:
    d = pd.read_table(f, sep = ' ', index_col = 0)
    dm = d.iloc[:,4:d.shape[1]].abs() # since coefficient sign doesn't mean anything with 0/1 SNPs, take abs. value
    dm = dm.mean(axis=0)
    print('Min: ' + str(dm.min()))
    print('Max: ' + str(dm.max()))
    print('Mean: ' + str(dm.mean()))
    print('std: ' + str(dm.std()))
    dm = (dm-dm.min())/(dm.max()-dm.min())
    print(dm.head())

  elif "_imp" in f:
    dm = pd.read_table(f, sep='\t',header=None, index_col=0) 
    print('Min: ' + str(dm.min()))
    print('Max: ' + str(dm.max()))
    print('Mean: ' + str(dm.mean()))
    print('std: ' + str(dm.std()))
    dm = (dm-dm.min())/(dm.max()-dm.min())
    print(dm.head())

  else:
    print("Predicted value file not recognized: %s" % f)
  
  if count == 1:
    merged = pd.DataFrame(dm)
  else:
    #d = pd.concat([d, pd.DataFrame(j)], axis=1)
    merged = pd.merge(merged, pd.DataFrame(dm), how='inner', left_index=True, right_index=True)
  count += 1

print('\n\nEnsemble:')
ensemble = merged.mean(axis=1)
print(ensemble.head())



ensemble.to_csv(SAVE, sep='\t', header=False, index=True)
