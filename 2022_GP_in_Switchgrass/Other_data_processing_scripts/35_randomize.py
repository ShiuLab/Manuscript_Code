import os,sys
import pandas as pd
import numpy
from numpy import random
from numpy.random import shuffle
file = sys.argv[1]
df = pd.read_csv(file, sep=',', index_col = None, header = 0)
# rowname = df.index.tolist()
# row = shuffle(rowname)
# df.index = rowname
shuffle(df.ID)
df.to_csv(file, index=False, header=True,sep=",")
