###########################################################
# Sum the SHAP feature interaction values across instances
# 
#
# Written by: Peipei Wang
###########################################################

import sys,os,argparse
import datatable as dt
import pandas as pd
import numpy as np # Added by Kenia 02/09/2024
from glob import glob
from tqdm import tqdm
from scipy import sparse

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def main():
    parser = argparse.ArgumentParser(
        description='Sum the SHAP feature interaction values across instances.')

	### Input arguments ###
	# Required
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-path', help='path saving the interaction score files', required=True)
    req_group.add_argument('-save', help='name to save', required=True)
    # Optional (Added by Kenia 02/01/2023)
    opt_group = parser.add_argument_group(title='OPTIONAL INPUT')
    opt_group.add_argument('-y', help='name of label', default='')
    opt_group.add_argument('-dtype', help='name of feature type', default='')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    if not args.path.endswith('/'):
        path = args.path + '/'
    else:
        path = args.path
        
    n = 0
    if (args.y == '' and args.dtype == ''): # Modified by Kenia 02/01/2023
        # Get row and column feature names (.npz files don't have such info)
        for file in glob(args.path + 'shap_interaction_scores_*.txt'):
            idx = dt.fread(file).to_pandas()
            idx = idx.iloc[:,0]
            idx.name = 'ID'
        for file in tqdm(os.listdir(args.path)):
            if file.endswith('.npz') and file.startswith('shap_interaction_scores'):
                if n == 0:
                    # shap = pd.read_csv(args.path + file,header=0,index_col=0,sep='\t')
                    # SHAP = pd.DataFrame(0,index = shap.index,columns=shap.columns)
                    shap = sparse.load_npz(args.path + file)
                    shap = pd.DataFrame(shap.toarray(), index=idx, columns=idx)
                    SHAP = pd.DataFrame(0, index=shap.index, columns=shap.columns)
                    n = n + 1
                # shap = pd.read_csv(args.path + file,header=0,index_col=0,sep='\t')
                shap = sparse.load_npz(args.path + file)
                shap = pd.DataFrame(shap.toarray(), index=idx, columns=idx)
                SHAP = SHAP.add(shap.abs())
    else: 
        # Get row and column feature names (.npz files don't have such info)
        for file in glob(args.path + 'shap_interaction_scores_*.txt'):
            if (args.y in file and args.dtype in file):
                idx = dt.fread(file).to_pandas()
                idx = idx.iloc[:,0]
                idx.name = 'ID'
        for file in tqdm(os.listdir(args.path)):
            if (file.endswith('.npz') and file.startswith('shap_interaction_scores')) \
                and (args.y in file and args.dtype in file):
                if n == 0:
                    shap = sparse.load_npz(args.path + file)
                    shap = pd.DataFrame(shap.toarray(), index=idx, columns=idx)
                    SHAP = pd.DataFrame(0,index = shap.index,columns=shap.columns)
                    n = n + 1
                shap = sparse.load_npz(args.path + file)
                shap = pd.DataFrame(shap.toarray(), index=idx, columns=idx)
                SHAP = SHAP.add(shap.abs())

    # Commented out by Kenia 02/09/2024
    # Res = pd.DataFrame(columns = ['Feature1','Feature2','Interaction'])
    # for i in range(0,SHAP.shape[0]):
    #     for j in range(i+1,SHAP.shape[0]):
    #         Res.loc[Res.shape[0],:] = [SHAP.index[i],SHAP.index[j],SHAP.iloc[i,j]]
    #     print(i)
    # Added by Kenia 02/09/2024
    try: # Added by Kenia 05/03/2024
        Res = SHAP.where(np.triu(SHAP, k=1).astype(bool)).stack()
        Res.rename('Interaction', inplace=True)
        Res.rename_axis(index=['Feature1', 'Feature2'], inplace=True)
        Res.to_csv(args.save + '_summed.txt',sep='\t',header=True,index=True) # Modified by Kenia 02/09/2024
    except UnboundLocalError as e:
        print(e)
        print("Do not pass -y argument and ensure -path ends in / to avoid this error.")
        sys.exit(1)
    # Res2 = Res[Res['Interaction'] > 0] # Commented out by Kenia 02/01/2024
    # Res2.to_csv(args.save + '_non_zero.txt',sep='\t',header=True,index=True)

if __name__ == '__main__':
    main()
