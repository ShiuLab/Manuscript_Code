#!/usr/bin/env python3
"""
Table S8
Combine GO/pathway enrichment results into an excel file
Arguments:
    path: path to csv files
    pattern: csv file name pattern to match
    dtype: data type (SNP, ORF, CNV)
    atype: annotation type (go, pwy)
    save: path to save excel file to
Usage:
    python combine_go_pwy_enrichment.py <path> <pattern> <dtype> <atype>
Examples:
    python Table_S8_combine_go_pwy_enrichment.py \
        Scripts/Data_Vis/Section_3/PWY_Enrichment/SNPs_fs/ \
        ORA_PWY_Genes_ SNP pwy Scripts/Data_Vis/Section_3/
    python Table_S8_combine_go_pwy_enrichment.py \
        Scripts/Data_Vis/Section_3/PWY_Enrichment/PAVs_fs/ \
        ORA_PWY_Genes_[A-Z0-9]+_pav PAV pwy Scripts/Data_Vis/Section_3/
    python Table_S8_combine_go_pwy_enrichment.py \
        Scripts/Data_Vis/Section_3/PWY_Enrichment/CNVs_fs/ \
        ORA_PWY_Genes_[A-Z0-9]+_cnv CNV pwy Scripts/Data_Vis/Section_3/

    python Table_S8_combine_go_pwy_enrichment.py \
        Scripts/Data_Vis/Section_3/GO_Enrichment/SNPs_fs/ \
        ORA_Genes_ SNP go Scripts/Data_Vis/Section_3/
    python Table_S8_combine_go_pwy_enrichment.py \
        Scripts/Data_Vis/Section_3/GO_Enrichment/PAVs_fs/ \
        ORA_Genes_[A-Z0-9]+_pav PAV go Scripts/Data_Vis/Section_3/
    python Table_S8_combine_go_pwy_enrichment.py \
        Scripts/Data_Vis/Section_3/GO_Enrichment/CNVs_fs/ \
        ORA_Genes_[A-Z0-9]+_cnv CNV go Scripts/Data_Vis/Section_3/
Returns:
    excel file with all GO/pathway enrichment results
    e.g. SNP_pathway_enrichment.xlsx or SNP_GO_enrichment.xlsx
"""

import os
import sys
import re
import pandas as pd
from pathlib import Path

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Arguments
path = sys.argv[1]  # path to csv files
pattern = sys.argv[2]  # csv file name pattern to match
dtype = sys.argv[3]  # data type (SNP, ORF, CNV)
atype = sys.argv[4]  # annotation type (go, pwy)
save = sys.argv[5]  # path to save excel file to

# Pathway names
pwy_info = pd.read_csv(
    "../Co-function/Data/MetaCyc/All_pathways_S288c_names.txt", sep="\t")
print(pwy_info.columns)

# Create and save excel file
if (atype == "pwy" or atype == "pathway"):
    excel_file = Path(save) / f"{dtype}_pathway_enrichment.xlsx"
    writer = pd.ExcelWriter(excel_file)
    pd.DataFrame({}).to_excel(writer, sheet_name="All")  # create empty sheet
    writer.close()
    with pd.ExcelWriter(excel_file, mode="a", if_sheet_exists="overlay") as writer:
        for csvfilename in os.listdir(path):
            if re.search(pattern, str(csvfilename)):
                print(csvfilename)
                df = pd.read_csv(Path(path) / csvfilename, sep="\t")
                print(df.shape)
                # remove non-significant pathways
                df = df.loc[df["p.val"] <= 0.05, :]
                print(df.shape)
                df['PWY'] = df['PWY'].str.strip()
                df = df.merge(pwy_info, left_on="PWY", right_on="Object ID")
                sheet_name = os.path.splitext(csvfilename)[0].split("_")[3]
                df.to_excel(writer, sheet_name=sheet_name)
                df.insert(0, "Environment", sheet_name)
                df.to_excel(writer, sheet_name="All", index=False,
                            header=False, startrow=writer.sheets['All'].max_row)

if (atype == "go" or atype == "GO"):
    excel_file = Path(save) / f"{dtype}_GO_enrichment.xlsx"
    writer = pd.ExcelWriter(excel_file)
    pd.DataFrame({}).to_excel(writer, sheet_name="All")  # create empty sheet
    writer.close()
    with pd.ExcelWriter(excel_file, mode="a", if_sheet_exists="overlay") as writer:
        for csvfilename in os.listdir(path):
            if re.search(pattern, str(csvfilename)):
                print(csvfilename)
                df = pd.read_csv(Path(path) / csvfilename, sep="\t")
                print(df.shape)
                # remove non-significant pathways
                df = df.loc[df["p.val"] <= 0.05, :]
                print(df.shape)
                sheet_name = os.path.splitext(csvfilename)[0].split("_")[2]
                df.to_excel(writer, sheet_name=sheet_name)
                df.insert(0, "Environment", sheet_name)
                df.to_excel(writer, sheet_name="All", index=False,
                            header=False, startrow=writer.sheets['All'].max_row)
