##############################################################################
# Convert more common site to 1 and less common site to -1 for the SNP matrix, 
# and transpose the matrix
# 
# Written by: Peipei Wang
##############################################################################

import pandas as pd
import sys,os,argparse

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def main():
    parser = argparse.ArgumentParser(
        description='Convert more common site to 1 and less common site to -1 for the SNP matrix and transpose the matrix.')

	### Input arguments ###
	# Required
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-file', help='the SNP matrix file to be converted and transposed', required=True)
    req_group.add_argument('-save', help='name to save the output file', required=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    df_MAF = pd.read_csv(args.file,header=0,index_col=0,sep=',')
    remain = df_MAF.sum(axis=1)
    index = remain.index.tolist()
    df_convert = df_MAF.copy(deep=True)
    for i in range(0,remain.shape[0]):
        if remain.iloc[i] <= df_MAF.shape[1]*0.5: # more common is 0, less common is 1
            df_convert.loc[index[i],df_convert.loc[index[i],:]==1] = -1
            df_convert.loc[index[i],df_convert.loc[index[i],:]==0] = 1
        else: # more common is 1, less common is 0
            df_convert.loc[index[i],df_convert.loc[index[i],:]==0] = -1
    df_convert = df_convert.T
    df_convert.to_csv(args.save + '_converted.csv',index=True, header=True,sep=",")

if __name__ == '__main__':
    main()

