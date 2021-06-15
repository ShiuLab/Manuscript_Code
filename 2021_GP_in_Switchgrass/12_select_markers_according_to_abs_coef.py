import sys,os, argparse 
import pandas as pd
import numpy as np
from scipy.stats import zscore
def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for selecting top n markers based on the abs coef')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-coef', help='coef of markers', required=True)
	req_group.add_argument('-start', help='start number', required=True)
	req_group.add_argument('-stop', help='end number', required=True)
	req_group.add_argument('-step', help='step', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	coef_file = args.coef # coef file
	start_number = int(args.start)
	stop_number = int(args.stop)
	step = int(args.step)
	coef = pd.read_csv(coef_file,header=0,index_col=None,sep=',')
	res = pd.DataFrame(coef.median(axis=0))
	res['abs_coef'] = abs(res[0])
	res = res.sort_values(['abs_coef'], ascending=False)
	for i in range(start_number,stop_number,step):
		with open('Markers_top%s.txt'%(i), 'w') as filehandle:
			for m in res.index.tolist()[0:i]:
				filehandle.write('%s\n' % m)


if __name__ == '__main__':
	main()

