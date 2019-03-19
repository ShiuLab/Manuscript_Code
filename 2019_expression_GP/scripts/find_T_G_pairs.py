"""
Find SNP with greatest absolute value coefficient/importance score within the LD limit of a transcript 
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
  if sys.argv[i].lower() == "-ld":
    LD = int(sys.argv[i+1])
  if sys.argv[i].lower() == "-m":
    M = sys.argv[i+1]

print('Reading in coefficient scores...')
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

print('Processing transcript and SNP scores...')
# For results from rrBLUP & BGLR (i.e. output into accuracy, *_coef.csv, and *_yhat.csv)
if TYP.lower() == 'coef':
  # Get mean values for SNPs with location info available (i.e. starts with chr)
  #cols = [c for c in g.columns if c.lower()[:3] == 'chr'] # For SNP data with SNPs named ChrX_###
  cols = [c for c in g.columns if c[:1] == 'S'] # For SNP data with SNPs named SX_###
  g = g[cols]
  g_m = pd.DataFrame(g.mean())
  g_m.columns = ['coef']
  g_m['coef'] = g_m['coef'].abs()
  g_m['Full_Loc'] = g_m.index

  # Get mean values for transcripts with one-to-one mapping
  t_m = pd.DataFrame(t.mean())
  t_m.columns = ['coef']
  t_m['coef'] = t_m['coef'].abs()
  t_m = t_m[t_m.index.isin(transcripts_to_keep)]
  t_m['v3_gene_model'] = t_m.index.tolist()

# For *_imp results (i.e. from Machine-Learning Pipeline):
elif TYP.lower() == 'imp':
  # Get transcripts with one-to-one mapping
  t_m = t[t['v3_gene_model'].isin(transcripts_to_keep)]

  # Get SNPs with location available
  g_m = g[g['Full_Loc'].str.startswith('Chr')]
  print(g_m.head())
else:
  print('Need to specify data input type as coef or imp')

# Pull location info out of SNP name

g_m['chromo'], g_m['loc'] = g_m['Full_Loc'].str.split('_', 1).str
g_m['loc'] = pd.to_numeric(g_m['loc'])
g_m['chromo'] = g_m['chromo'].str.replace('S','Chr')
print(g_m.head(3))

# Merge with key and calculate gene center location
t_m = t_m.merge(k_dups, how='inner', on='v3_gene_model')
t_m['trans_mid'] = t_m[['v4_start', 'v4_end']].mean(axis=1).round()
t_m = t_m[~t_m['v4_chr'].astype(str).str.startswith('B73V4')]
print(t_m.head(3))

# For each chromosome find SNP in LD with greatest coef
print('Finding top SNP within %i of each gene...' % LD)
count = 0
r = pd.DataFrame()

for row in t_m.itertuples():
  count += 1
  loc_e = row.trans_mid + LD
  loc_s = row.trans_mid - LD
  loc_chromo = row.v4_chr
  SNPs_in_region = g_m[g_m.eval('(chromo == @loc_chromo) & (loc >= @loc_s) & (loc <= @loc_e)')]

  if M.lower() == 'top':
    try:
      top_SNP = SNPs_in_region['coef'].idxmax()
      top_SNP_loc = SNPs_in_region['loc'][top_SNP]
      top_SNP_coef = SNPs_in_region['coef'][top_SNP]
      r = r.append({'Chromo': row.v4_chr,'v3_gene_model': row.v3_gene_model, 'gene_s': row.v4_start, 'gene_e':row.v4_end,
        'gene_mid':row.trans_mid, 'gene_coef':row.coef, 'SNP': top_SNP, 'SNP_loc': top_SNP_loc, 
        'SNP_coef': top_SNP_coef}, ignore_index=True)
    except:
      r = r.append({'Chromo': row.v4_chr, 'v3_gene_model': row.v3_gene_model, 'gene_s': row.v4_start, 'gene_e':row.v4_end,
        'gene_mid':row.trans_mid, 'gene_coef':row.coef, 'SNP': 'na', 'SNP_loc': 'na', 
        'SNP_coef': 'na'}, ignore_index=True)
  
  elif M == '95th':
    try:
      top_SNP = SNPs_in_region['coef'].idxmax()
      top_SNP_loc = SNPs_in_region['loc'][top_SNP]
      x95th_SNP_coef = SNPs_in_region['coef'].quantile(0.95)
      r = r.append({'Chromo': row.v4_chr,'v3_gene_model': row.v3_gene_model, 'gene_s': row.v4_start, 'gene_e':row.v4_end,
        'gene_mid':row.trans_mid, 'gene_coef':row.coef, 'SNP': top_SNP, 'SNP_loc': top_SNP_loc, 
        'SNP_coef': x95th_SNP_coef}, ignore_index=True)
    except:
      r = r.append({'Chromo': row.v4_chr, 'v3_gene_model': row.v3_gene_model, 'gene_s': row.v4_start, 'gene_e':row.v4_end,
        'gene_mid':row.trans_mid, 'gene_coef':row.coef, 'SNP': 'na', 'SNP_loc': 'na', 
        'SNP_coef': 'na'}, ignore_index=True)
  if count % 1000 == 0:
    print('\tFinished %i genes' % count)

print(r.head())

SAVE = SAVE + ".csv"
r.to_csv(SAVE, sep=',')
