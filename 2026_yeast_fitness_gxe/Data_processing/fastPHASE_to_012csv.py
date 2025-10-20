"""This script converts the output of the fastPHASE_to_csv.py script, where
SNP genotypes are encoded as -1 (homozygous major), 0 (heterozygous), and 1
(homozygous minor), into a CSV file where genotypes are encoded as 0 (homozygous ref),
1 (heterozygous), and 2 (homozygous alt).

Example Usage:
	path=/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/0_raw_data
	file=1011Matrix_diploid750_biallelic_maf0.05_09102025_genotypes.csv
	frq=1011Matrix_diploid750_biallelic_maf0.05_09102025.fastphase.frq
	vcf=1011Matrix_diploid750_biallelic_maf0.05_11032022.recode.vcf
	save=geno_012.csv
	python fastPHASE_to_012csv.py -path $path -file $file -frq $frq -vcf $vcf -save $save


Args:
	path (str): Path to working directory.
	file (str): Name of fastPHASE_to_csv.py output file in the working directory.
	frq (str): Name of PLINK .frq file containing which allele is major/min.
	vcf (str): Name of VCF file containing the ref/alt allele genotypes.
	save (str): Name to save genotype CSV file as.
"""

import os
import argparse
import datatable as dt
import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert -1/0/1 genotypes based on major/minor alleles to 0/1/2 genotype CSV format."
    )
    parser.add_argument("-path", type=str, required=True,
                        help="Path to working directory.")
    parser.add_argument("-file", type=str, required=True,
                        help="Name of fastPHASE_to_csv.py output file in the working directory.")
    parser.add_argument("-frq", type=str, required=True,
                        help="Name of PLINK .frq file containing which allele is major/min.")
    parser.add_argument("-vcf", type=str, required=True,
                        help="Name of VCF file containing the ref/alt allele genotypes.")
    parser.add_argument("-save", type=str, required=True,
                        help="Name to save genotype CSV file as.")
    return parser.parse_args()


def convert_to_012(orig_geno, frq_dict, vcf_dict, save=""):
    """ Convert -1/0/1 genotypes to 0/1/2 genotypes based on major/minor and
    ref/alt allele nucleotide information."""

    # Create an empty dictionary to store the converted genotypes
    converted_geno = {}

    # Iterate through each SNP column
    for snp in tqdm(orig_geno.columns, desc="Converting genotypes to 0/1/2 format"):
        # Fetch the major/minor and ref/alt alleles for the SNP
        minor_allele = frq_dict[snp]["minor"]
        major_allele = frq_dict[snp]["major"]
        ref_allele = vcf_dict[snp]["REF"]
        alt_allele = vcf_dict[snp]["ALT"]

        # Determine the mapping from major/minor to ref/alt
        # -1: homozygous major, 0: heterozygous, 1: homozygous minor
        if (major_allele == ref_allele) and (minor_allele == alt_allele):
            genotype_mapping = {-1: 0, 0: 1, 1: 2}  # major=ref, minor=alt
        elif (major_allele == alt_allele) and (minor_allele == ref_allele):
            genotype_mapping = {-1: 2, 0: 1, 1: 0}  # major=alt, minor=ref

        # Apply the mapping to the genotype column
        converted_geno[snp] = orig_geno[snp].map(genotype_mapping)

    # Convert the dictionary back to a DataFrame
    converted_geno = pd.DataFrame(converted_geno, index=orig_geno.index)
    converted_geno.index.name = "ID"

    print(f"Saving converted genotype data to {save}...")
    converted_geno.reset_index(inplace=True)
    converted_geno = dt.Frame(converted_geno)
    converted_geno.to_csv(save, quoting="none")

    return converted_geno


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    path = args.path
    file = args.file
    frq = args.frq
    vcf = args.vcf
    save = args.save

    # Change to working directory
    os.chdir(path)

    # Load the input files
    orig_geno = dt.fread(file).to_pandas()  # -1,0,1 genotypes (major/minor)
    frq_data = dt.fread(frq).to_pandas()
    vcf_data = dt.fread(vcf, skip_to_line=59).to_pandas()

    # For faster lookup, create dictionaries of major and minor alleles
    # example format of dictionary: {'chr_pos': {'minor': 'A', 'major': 'G'}}
    frq_data["ID"] = frq_data["SNP"].str.replace(":", "_")  # change to chr_pos
    frq_data.rename(columns={"A1": "minor", "A2": "major"}, inplace=True)
    frq_dict = frq_data.set_index(
        "ID")[["minor", "major"]].to_dict(orient='index')

    # Create a dictionary of reference and alternate alleles from the VCF
    # example format of dictionary: {'chr_pos': {'ref': 'G', 'alt': 'A'}}
    vcf_data.insert(0, "SNP_ID", vcf_data["#CHROM"].astype(str).str.strip() +
                    "_" + vcf_data["POS"].astype(str).str.strip())
    vcf_dict = vcf_data.set_index("SNP_ID")[["REF", "ALT"]].\
        to_dict(orient='index')

    # Now, convert the -1/0/1 genotypes to 0/1/2 genotypes
    orig_geno.set_index("index", inplace=True)  # set sample IDs as index
    converted_geno = convert_to_012(orig_geno, frq_dict, vcf_dict, save=save)
