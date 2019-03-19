"""
Find cis and trans eQTL for each transcript with the greatest absolute value coefficient/importance score 
"""
import sys, os
import pandas as pd

key_file = '/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/00_RawData/v3_v4_xref.txt'
g_coef_file = '/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/03_rrBLUP/test_geno_coef.csv'
t_coef_file = '/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/03_rrBLUP/test_trans_coef.csv'
LD = 1000
TYP = 'coef' 
SAVE = 'test_imp_compare'
M = 'top'

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-g":
    g_coef_file = sys.argv[i+1]
  if sys.argv[i] == "-t":
    t_coef_file = sys.argv[i+1]
  if sys.argv[i] == "-key":
    key_file = sys.argv[i+1]
  if sys.argv[i] == "-type":
    TYP = sys.argv[i+1]
  if sys.argv[i] == "-save":
    SAVE = sys.argv[i+1]
  if sys.argv[i].lower() == "-cis":
    CIS = sys.argv[i+1]
  if sys.argv[i].lower() == "-trans":
    TRANS = sys.argv[i+1]
  if sys.argv[i].lower() == "-m":
    M = sys.argv[i+1]

print('Proceessing cis- & trans-eQTL results')
cis = pd.read_csv(CIS, sep='\t', header=0, usecols=['SNP','gene','beta','t-stat','p-value','FDR'])
trans = pd.read_csv(TRANS, sep='\t', header=0, usecols=['SNP','gene','beta','t-stat','p-value','FDR'])
cis = cis[cis.FDR <= 0.05]
trans = trans[trans.FDR <= 0.05]
print('\nNumber of significant cis-eQTL %i' % len(cis.index))
print('Number of significant trans-eQTL %i' % len(trans.index))
eqtl = pd.concat([cis, trans])
print('Total number of significant eQTL %i' % len(eqtl.index))

print('\nRead and processing in coefficients/importance scores...')
#g = pd.read_csv(g_coef_file, sep=' ', header=0)
if TYP.lower() == 'coef':
  g = pd.read_csv(g_coef_file, sep=' ', header=0)
  t = pd.read_csv(t_coef_file, sep=' ', header=0)
elif TYP.lower() == 'imp':
  g = pd.read_csv(g_coef_file, sep='\t', header=None, names=['Full_Loc', 'coef'])
  t = pd.read_csv(t_coef_file, sep='\t', header=None, names=['v3_gene_model', 'coef'])
else:
  print('Need to specify data input type as coef or imp!!!')
  quit()

# Get list of transcripts with one-to-one mapping between AGP v3 and v4
k = pd.read_csv(key_file, sep='\t', header=0, usecols=['v3_gene_model', 'v4_gene_model', 'v4_chr', 'v4_start', 'v4_end'])
k_dups = k.drop_duplicates(subset=['v3_gene_model','v4_gene_model'], keep='first')
k_dups = k_dups.drop_duplicates(subset='v3_gene_model', keep=False)
k_dups = k_dups.drop_duplicates(subset='v4_gene_model', keep=False)
transcripts_to_keep = k_dups['v3_gene_model'].tolist()

# For results from rrBLUP & BGLR (i.e. output into accuracy, *_coef.csv, and *_yhat.csv)
if TYP.lower() == 'coef':
  # Get mean values for SNPs with location info available (i.e. starts with chr)
  #cols = [c for c in g.columns if c.lower()[:3] == 'chr'] # For SNP data with SNPs named ChrX_###
  #cols = [c for c in g.columns if c[:1] == 'S'] # For SNP data with SNPs named SX_###
  #g = g[cols]
  g_m = pd.DataFrame(g.mean())
  g_m.columns = ['coef']
  g_m['coef'] = g_m['coef'].abs()
  g_m['Full_Loc'] = g_m.index

  # Get mean values for transcripts with one-to-one mapping
  t_m = pd.DataFrame(t.mean())
  t_m.columns = ['coef']
  t_m['coef'] = t_m['coef'].abs()
  t_m['v3_gene_model'] = t_m.index.tolist()

# For *_imp results (i.e. from Machine-Learning Pipeline):
elif TYP.lower() == 'imp':
  # Get transcripts with one-to-one mapping
  t_m = t
  t_m.index = t_m['v3_gene_model']

  # Get SNPs with location available
  g_m = g
  g_m.index = g_m['Full_Loc']
else:
  print('Need to specify data input type as coef or imp')



g_m['Full_Loc'] = g_m['Full_Loc'].str.replace('B73V4_','')
g_m['chromo'], g_m['loc'] = g_m['Full_Loc'].str.split('_', 1).str

g_m['loc'] = pd.to_numeric(g_m['loc'])
g_m['chromo'] = g_m['chromo'].str.replace('S','')
g_m['chromo'] = g_m['chromo'].str.replace('Chr','')
g_m['chromo'] = g_m['chromo'].str.replace('ctg','')

print(g_m.head())
print(t_m.head())


# Pull top cis and trans eQTL for each transcript with a significant eQTL!

df_cols = columns=['v3_gene_model', 'gene_imp', 'top_eqtl', 'eqtl_imp']
eqtl_list = []

for e in eqtl.gene.unique():
  t_imp = t_m.loc[e]['coef']
  top_snp_imp = 0
  top_snp = ''
  eqtl_tmp = eqtl[eqtl.gene == e]
  all_eqtl_options = pd.DataFrame()
  
  for eqtl_x in eqtl_tmp['SNP']:
    try:
      eqtl_x_chromo, eqtl_loc = eqtl_x.strip().split('_')
    except:
      xx, eqtl_x_chromo, eqtl_loc = eqtl_x.strip().split('_')
    eqtl_x_chromo = ''.join([i for i in eqtl_x_chromo if i.isdigit()])
    loc_e = float(eqtl_loc) + LD
    loc_s = float(eqtl_loc) - LD
    all_eqtl_options_tmp = g_m[g_m.eval('(chromo == @eqtl_x_chromo) & (loc >= @loc_s) & (loc <= @loc_e)')]
    all_eqtl_options = pd.concat([all_eqtl_options, all_eqtl_options_tmp])
  
  all_eqtl_options = all_eqtl_options.sort_values(by='coef', ascending=False)
  top_snp = all_eqtl_options['Full_Loc'].iloc[0]
  top_snp_imp = all_eqtl_options['coef'].iloc[0]

  eqtl_list.append([e, t_imp, top_snp, top_snp_imp])


eqtl_df = pd.DataFrame(eqtl_list, columns=df_cols)
print('Snapshot of eQTL results:')
print(eqtl_df.head())



# Save results
eqtl_SAVE = SAVE + "_eQTL.csv"
eqtl_df.to_csv(eqtl_SAVE, sep=',')


print('Done!')