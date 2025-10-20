#!/usr/bin/env python3

"""
Convert fastPHASE v1.4.8 imputed biallelic genotypes file to CSV matrix format

Note: The minor allele of the fastPHASE file is encoded as 1 and 
    the major allele is encoded as 2.

Example:
    # path=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/
    path=/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/0_raw_data
    file=1011Matrix_diploid750_biallelic_maf0.05_11032022_genotypes.out
    # map=1011Matrix_diploid750_biallelic_maf0.05_11032022.plink.map
    map=1011Matrix_diploid750_biallelic_maf0.05_09102025.fastphase.frq
    recoded=1011Matrix_diploid750_biallelic_maf0.05_11032022.fastphase.chr-0.recode.phase.inp
    # save=1011Matrix_diploid750_biallelic_maf0.05_11032022_genotypes.csv
    save=1011Matrix_diploid750_biallelic_maf0.05_09102025_genotypes.csv
    # python fastPHASE_to_csv.py -path $path -file $file -map $map -save $save
    python fastPHASE_to_csv.py -path $path -file $file -map $map -recoded $recoded -save $save
    
Args:
    path (str): Path to working directory.
    file (str): Name of fastPHASE imputed _genotype.out file.
    map (str): Name of PLINK map file containing chromosome and bp position of genotypes.
    recoded (str): Name of the PLINK output file after recoding into major/minor alleles in the fastphase format.
    save (str): Name to save matrix file as.

Returns:
    Saves a comma-separate values file in the working directory. Genotypes are
    encoded as -1 (homozygous major), 0 (heterozygous), 1 (homozygous minor)[1].

Reference:
    [1] Kim, S., & Lee, C. Y. (2020). A consistent approach to the genotype 
        encoding problem in a genome-wide association study of continuous 
        phenotypes. PloS one, 15(7), e0236139.
"""

__author__ = "Kenia Segura Abá"

import sys
import os
import argparse
import warnings
import pandas as pd
import datatable as dt
from tqdm import tqdm
import pyspark


def warn(*args, **kwargs):
    pass


warnings.warn = warn


def fastPHASE_to_csv(inf, recoded, map, save):
    """ Convert fastPHASE imputed biallelic genotypes file to CSV """
    # read the PLINK recoded genotypes file to get the genotype order
    recoded_inl = recoded.readline()
    while not recoded_inl.strip().startswith("P"):
        recoded_inl = recoded.readline()

    # fetch the snp IDs based on the order in the recoded file
    recoded_inl = pd.Series(
        recoded_inl.strip().split(" ")[1:])  # skip the "P" at index 0
    recoded_inl.name = "P"
    assert recoded_inl.eq(map["PLINK_ID"]).all()  # sanity check, returns True!

    # read the imputed genotypes line by line
    inl = inf.readline()  # read the first line of the file

    geno = {}  # dict of IDs and recoded genotypes

    with tqdm(total=751) as pbar:  # progress bar
        while inl.strip() != "END GENOTYPES":
            if inl.strip() == "BEGIN GENOTYPES":
                pass

            if inl.strip().startswith("#"):
                # fetch individual IDs
                ID = inl.strip().split(" ")[2]  # individual ID
                if ID not in geno.keys():  # check if ID exists
                    geno[ID] = {}

                # fetch genotypes for all SNPs
                inl = inf.readline().strip().split(" ")  # allele 1
                inl2 = inf.readline().strip().split(" ")  # allele 2

                # recode SNP genotype
                for i in range(len(inl)):
                    # SNP name ---------------------------------This is where I went wrong
                    snp = map.iloc[i, 6]  # fetch SNP ID in chr_pos format
                    if snp not in geno[ID].keys():
                        if inl[i] == "2" and inl2[i] == "2":  # homozygous major
                            geno[ID][snp] = -1
                        elif inl[i] != inl2[i]:  # heterozygous
                            geno[ID][snp] = 0
                        elif inl[i] == "1" and inl2[i] == "1":  # homozygous minor
                            geno[ID][snp] = 1
            inl = inf.readline()
            pbar.update(1)
    pbar.close()

    # save recoded genotype matrix to file
    geno = pd.DataFrame(geno)
    geno = geno.T
    geno.rename(columns={'index': 'ID'}, inplace=True)
    geno.reset_index(inplace=True)
    geno = dt.Frame(geno)
    geno.to_csv(save, quoting="none")
    return (0)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(
        description="Test script to convert vcf to matrix format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-path", required=True, type=str,
                           help="Path to working directory.")
    req_group.add_argument("-file", required=True, type=str,
                           help="Name of fastPHASE imputed _genotype.out file in the working directory.")
    req_group.add_argument("-map", required=True, type=str,
                           help="Name of PLINK .frq file containing chromosome and bp position of genotypes. It contains the SNP ID order in the recoded file.")
    req_group.add_argument("-recoded", required=True, type=str,
                           help="Name of the PLINK output file after recoding into major/minor alleles in the fastphase format.")
    req_group.add_argument("-save", required=True, type=str,
                           help="Name to save matrix file as.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()

    # Set working directory
    os.chdir(args.path)

    # Read in the recoded fastPHASE file, imputed genotypes file, & PLINK mapping file
    inf = open(args.file, "r")
    map = dt.fread(args.map).to_pandas()
    recoded = open(args.recoded, "r")

    # Add the SNP IDs to the map file
    map["ID"] = map["SNP"].str.replace(":", "_")  # change to chr_pos format
    map["PLINK_ID"] = map["SNP"].str.split(":", expand=True)[1]

    # Convert fastPHASE imputed genotypes to CSV
    fastPHASE_to_csv(inf=inf, recoded=recoded, map=map,
                     save=args.save)  # convert to CSV
    inf.close()
