###########################################################
# Sum the SHAP feature interaction values across instances
# 
#
# Written by: Peipei Wang
###########################################################

import sys,os,argparse
import pandas as pd

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
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    if not args.path.endswith('/'):
        path = args.path + '/'
    else:
        path = args.path
        
    n = 0
    for file in os.listdir(args.path):
        if file.endswith('.txt') and file.startswith('shap_interaction_scores'):
            if n == 0:
                shap = pd.read_csv(args.path + file,header=0,index_col=0,sep='\t')
                SHAP = pd.DataFrame(0,index = shap.index,columns=shap.columns)
                n = n + 1
            shap = pd.read_csv(args.path + file,header=0,index_col=0,sep='\t')
            SHAP = SHAP.add(shap.abs())

    Res = pd.DataFrame(columns = ['Feature1','Feature2','Interaction'])
    for i in range(0,SHAP.shape[0]):
        for j in range(i+1,SHAP.shape[0]):
            Res.loc[Res.shape[0],:] = [SHAP.index[i],SHAP.index[j],SHAP.iloc[i,j]]
        print(i)

    Res = Res.sort_values(by='Interaction',ascending=False)
    Res.to_csv(args.save + '.txt',sep='\t',header=True,index=True)

    Res2 = Res[Res['Interaction'] > 0]
    Res2.to_csv(args.save + '_non_zero.txt',sep='\t',header=True,index=True)

if __name__ == '__main__':
    main()
