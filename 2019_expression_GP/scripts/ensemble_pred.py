"""
Take multiple _scores.txt and _yhat.csv files and 
generate an ensemble regression accuracy and yhats

Approach:
  1. For every replicate (default 100), grab the yhats from each model 
  2. Calculate the mean yhat for that replicate
  3. Calculate the correlation(ensemble yhat ~ y) and save in accuracy.txt
  4. Save all those ensemble yhats in the SAVE file

python ensemble_pred.py SAVE TRUE_Y_FILE Y_COL_NAME DTYPE_ID FILES...

After first 4 args, can give as many _scores.txt and _yhat.csv files as desired

"""


import sys, os
import pandas as pd
import numpy as np
import time

items = []

n = 100

SAVE = sys.argv[1]
PHENO = sys.argv[2]
TRAIT = sys.argv[3]
DTYPE = sys.argv[4]

# Read in scores from models to ensemble
for i in range (5,len(sys.argv)):
  items.append(sys.argv[i])

print(items)
# Read in true phenotypes and get trait column
pheno = pd.read_csv(PHENO, index_col=0)
pheno = pheno[TRAIT]


# For each replicate
count_x = 1
for x in range(1,n+1):
  start_time = time.time()
  # Calculate or pull mean yhat from each model
  count = 1
  yhat_name = 'cv_' + str(x)
  scores_name = 'rep_' + str(x)
  
  for f in items:
    if "_yhat.csv" in f:
      d = pd.read_csv(f, index_col = 0)
      d = d.drop_duplicates()
      dm = d[d.index.str.endswith(yhat_name)]
      dm = dm.transpose()

    elif "_scores.txt" in f:
      d = pd.read_table(f, sep='\t',header=0, index_col=0)
      dm = d[scores_name]

    else:
      print("Predicted value file not recognized: %s" % f)
    
    if count == 1:
      d = pd.DataFrame(dm)
    else:
      #d = pd.concat([d, pd.DataFrame(j)], axis=1)
      d = pd.merge(d, pd.DataFrame(dm), how='inner', left_index=True, right_index=True)
    count += 1
  ensemble = d.mean(axis=1)
  

  if count_x == 1:
    d_final = pd.DataFrame(ensemble)
  else:
    d_final = pd.merge(d_final, pd.DataFrame(ensemble), how='inner', left_index=True, right_index=True)

  ensemble = pd.merge(left = pd.DataFrame(ensemble), right = pd.DataFrame(pheno), left_index=True, right_index=True)
  cor = np.corrcoef(ensemble[0], ensemble[TRAIT])

  total_time = (time.time() - start_time)
  
  with open('accuracy.txt', 'a') as out:
    out.write('ensem\t%s\t%s\t%s\t%f\t%s\n' % (DTYPE, TRAIT, yhat_name, cor[0,1], total_time))
  count_x += 1



print(d_final.head())
d_final.to_csv(SAVE, header=True, index=True)
