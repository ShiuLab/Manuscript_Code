"""
Convert mapping results into ML like dataframe
"""

import pandas as pd
import numpy as np
import sys

for i in range (1,len(sys.argv),2):
    if sys.argv[i] == "-mapping":
    	m = sys.argv[i+1]

# Set up dataframe
mapped = pd.read_csv(m, sep='\t', index_col = 0)
index = np.unique(mapped.index)
columns = np.unique(mapped.Motif)
df = pd.DataFrame(index=index, columns=columns)

# Replace Nan with a 1 for every hit in the file
for l in open(m, 'r').readlines():
  if l.startswith("Sequence"):
    pass
  else:
    seq, motif = l.strip().split("\t")[0:2]
    df.ix[seq,motif] = 1

# Fix df headers & row names
df = df.fillna(0)
df.index = df.index.map(lambda x: str(x)[:-3])
df.columns = df.columns.map(lambda x: str(x)[:-8])

name = m + "_df.txt"
df.to_csv(name, sep = "\t")