#!/bin/python
"""
Match SNPs (binary, MAF > 0.05) to Arabidopsis thaliana TAIR10 reference genes
Written by: Kenia Segura Abá

Input:
    genotype file
    gff file
"""

import sys,os,argparse
import datatable
import pandas as pd
from alive_progress import alive_bar

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Match SNPs to Arabidopsis thaliana TAIR10 reference genes")
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument("-geno", help="genotype csv file", required=True)
    req_group.add_argument("-gff", help="reference genome gff file", required=True)
    req_group.add_argument("-save", help="save name for output file", required=True)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Read input files
    print("Reading in data...")
    geno = datatable.fread(args.geno).to_pandas()
    gff = open(args.gff, "r").readlines()

    # List of bi-allelic snps from diploid yeast isolates
    snps = geno.columns[1:] # exclude "ID" header

    # Extract all the genes and their information from the reference genome
    G = {} # directory to hold gene names, start, stop, and strand
    print("Extracting reference genes...")
    with alive_bar(len(gff), bar = "circles", spinner = "dots_waves") as bar: # progress bar
        for inl in gff: # loop through each line in gff file
            if inl.startswith("#"): # skip these lines
                pass
            elif inl.startswith(">"): # exit loop when genome sequence reached
                break
            else:
                tem = inl.split('\t') # tab delimeted elements are split into a list
                chr = tem[0] # the first element is the chromosome number location of the sequence
                type = tem[2] # the third element is the sequence type
                if chr not in G: # dictionary for each chromosome number
                    G[chr] = {}
                if type == "gene": # search genes
                    name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1].strip() # gene name
                    start = tem[3] # gene start position
                    stop = tem[4] # gene stop position
                    strand = tem[6] # forward (+) or reverse (-)
                    G[chr][name] = [start, stop, strand]
            bar()
 
    out = open("all_genes_TAIR10.txt", "w") # write TAIR10 genes to file
    for key, value in G.items(): # key is chromosome
        for k, v in value.items(): # k is gene name
            out.write("%s,%s,%s,%s,%s\n" % (key, k, v[0], v[1], v[2])) # v is start, stop, strand
    out.close()

    # Match SNPs to genes
    map = {}
    print("Mapping snps to reference genes...")
    with alive_bar(len(snps), bar = "circles", spinner = "dots_waves") as bar: # progress bar
        for s in snps: # loop through bi-allelic snps
            chr = s.split("_")[0] # chromosome number
            pos = s.split("_")[1] # position of snp
            for CHR, value in G.items(): # loop through reference genes
                for gene, v in value.items():
                    start = v[0]
                    stop = v[1]
                    if (chr == "Chr1" and CHR == "Chr1"): # chromosomes match
                        if (pos >= start and pos <= stop): # snp falls within genic region
                            map[s] = [gene] # map gene to snp
                        else:
                            pass
                    elif (chr == "Chr2" and CHR == "Chr2"):
                        if (pos >= start and pos <=stop):
                            map[s] = [gene]
                        else:
                            pass
                    elif (chr == "Chr3" and CHR == "Chr3"):
                        if (pos >= start and pos <=stop):
                            map[s] = [gene]
                        else:
                            pass
                    elif (chr == "Chr4" and CHR == "Chr4"):
                        if (pos >= start and pos <=stop):
                            map[s] = [gene]
                        else:
                            pass
                    elif (chr == "Chr5" and CHR == "Chr5"):
                        if (pos >= start and pos <=stop):
                            map[s] = [gene]
                        else:
                            pass
            bar()

    print("Saving file...")
    map = pd.DataFrame(map, index=[0])
    map = map.transpose()
    map.to_csv(args.save)

if __name__ == '__main__':
    main()
