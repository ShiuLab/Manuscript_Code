import sys, os, argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
	description="Convert feature file(s) into numpy data files to use as input to DNN models",
	epilog='https://github.com/azodichr')

parser.add_argument('-input', help='ML input file to convert', required=True)
parser.add_argument('-method', help='How to convert (options: omic_stack, ', required=True)
parser.add_argument('-y_name', help='Name of column you want to predict', default='pass')
parser.add_argument('-out', help='Save name for output numpy dataframe', default='default')
parser.add_argument('-sep', help='Seperator in input file', default='\t')

a = parser.parse_args()

input_df = pd.read_csv(a.input, index_col=0, sep=a.sep)

if a.out == 'default':
	a.out = a.input + "_" + 'a.method'

if a.y_name != 'pass':
	y = input_df[a.y_name].values
	np.save(a.out + '_y.npy', y)
	print('Final shape of y: %s' % str(y.shape))
	input_df.drop(labels=a.y_name, axis='columns', inplace=True)

if a.method.lower() == 'omic_stack':
	print('Stacking based on omic type...')

	features = input_df.columns.tolist()
	omic_data = []
	motifs = []
	for f in features:
		f_items = f.strip().split('_')
		if len(f_items) == 2:
			omic_data.append(f_items[0])
			motifs.append(f_items[1])
		else:
			input_df.rename(columns={f: 'PA_' + f}, inplace=True)
			omic_data.append('PA')
			motifs.append(f)
	omic_data = set(omic_data)
	motifs = set(motifs)

	print('Each instance will have %i columns (i.e. motifs) and %i rows (i.e. omic data types)' % (len(motifs), len(omic_data)))
	print('\n## Note the rows will be sorted alphabetically as shown below##\n')
	# Convert to multi-index df
	multi_index_list = input_df.columns.str.split('_', expand=True).values
	input_df.columns = pd.MultiIndex.from_tuples([(x[1], x[0]) for x in multi_index_list])
	input_df = input_df.sort_index(axis=1)

	# Convert to 3D numpy array
	# m = motifs, v = variations, n = instances
	count = 0
	for feat in list(input_df.columns.levels[0]):
		if count == 0:
			print(input_df[feat].head(5))
			out_array = input_df[feat].values
		out_array = np.dstack((out_array, input_df[feat].values))
		count += 1

print('\nSnapshot of first motfis on top 10 genes')
print(out_array[0:10, :, 0])
print('\nFinal shape of x: %s' % str(out_array.shape))
np.save(a.out + '_x.npy', out_array)

print('\n\nDone!')