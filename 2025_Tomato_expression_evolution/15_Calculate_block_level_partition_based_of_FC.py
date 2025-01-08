#this is a code modified from calculate_block_level_partitions.py after discussion with Shinhan
#calculate block level partition based of Fold change values

import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from matplotlib.colors import ListedColormap
import seaborn as sns
from pandas import DataFrame, concat
#_______________________________________________________________________________________________________________________#
#functions
def add_row(df, na_rows):
    for i in na_rows.index:
        line = DataFrame(np.nan, index=[f"no_synt_{i}"], columns=df.columns)
        df = concat([df.iloc[:i], line, df.iloc[i:]])
    return df
#_______________________________________________________________________________________________________________________#
#getting the alignmnts with median KS values between 0.5871012 and 0.89355>> evolved from the most recent GT event
ks_df = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/meadian_ks_for_syntanic_block.txt")
sub_ks_df = ks_df[ks_df["median_KS"] > 0.5871012]
sub_ks_df =  sub_ks_df[sub_ks_df["median_KS"] < 0.89355] 
algnment_list = list(sub_ks_df["block"])
tissue_cluster_groups = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/Tissue_cluster_groups_TPM.txt",header=0)
##_______________________________________________________________________________________________________________________#
ks_file = "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/ks_ka_tomato_paralogs_MCScanx_default"
tadem_file = "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/ITAG4.0_Tomota.tandem"
TPM_df = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/TPM_table_non_T_tissues.txt",header=0, index_col=0)
#remove Mature.anthers from the columns
TPM_df = TPM_df.drop("Mature.anthers", axis=1)

tissue_cluster_groups = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/Tissue_cluster_groups_TPM.txt",header=0)
#keep rows where tissue is equal to column names of TPM_df
tissue_cluster_groups = tissue_cluster_groups[tissue_cluster_groups["tissue"].isin(TPM_df.columns)]
#drop mature anthers from the tissue_cluster_groups
tissue_cluster_groups = tissue_cluster_groups[tissue_cluster_groups["tissue"] != "Mature.anthers"]
#change seedling to Root
tissue_cluster_groups["group"] = tissue_cluster_groups["group"].replace("Seedling","Root")
#change the order of group column to "Fruit","Seed","Root","Callus","Floral","Stem","Trichome" order
group_order = tissue_cluster_groups["group"].to_list()

custom_order = {'Fruit': 1, 'Seed': 2, 'Root': 3, 'Callus': 4, 'Floral': 5, 'Stem': 6, 'Trichome': 7}

sorted_group_order = sorted(group_order, key=lambda x: custom_order.get(x, float('inf')))

#order rows of tissue_cluster_groups based on the custom_order
tissue_cluster_groups["group"] = tissue_cluster_groups["group"].map(custom_order)
#sort the rows based on the group column
tissue_cluster_groups = tissue_cluster_groups.sort_values("group")

#change the order of columns in TPM_df to match the order of rows in tissue_cluster_groups
TPM_df = TPM_df[tissue_cluster_groups["tissue"]]
                                                                                                                                                          
tadem_df =  pd.read_csv(tadem_file, header=None)
annotation_file = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/extracted_gene_annotation.txt", header=0, index_col=0)
#out = open("meadian_ks_for_syntanic_block.txt","w")
     
#_______________________________________________________________________________________________________________________#
Block_dict = {}
in_list = open(ks_file,"r").readlines()[11:]
for lines in in_list:      
    if lines.startswith("## Alignment"):
        almn_nm = lines.split(":")[0].replace("## ","").lstrip().rstrip()
        #val_list = []
        #print(almn_nm)
        if almn_nm in algnment_list:
            if almn_nm not in Block_dict.keys():
                Block_dict[almn_nm] = {}
                Block_dict[almn_nm]["Left"] = []
                Block_dict[almn_nm]["Right"] = []
                
                
            print("alignmnt found")
            switch = "yes"
        else:
            switch = "no"
    else:
        if switch == "yes":
            left_gene = lines.strip("").split("\t")[1].lstrip().rstrip()
            right_gene = lines.strip("").split("\t")[2].lstrip().rstrip()
            Block_dict[almn_nm]["Left"].append(left_gene)
            Block_dict[almn_nm]["Right"].append(right_gene)
(Block_dict.keys())
 #_______________________________________________________________________________________________________________________#
#trshould to call two gens are differentially expressed based on log2(TPM_gene1/TPM_gene2) values
diff_thre = 10

#dict to save partition scores
partition_score_dict = {}
#dict to save partition direction (this is numbers of genes partitioned in perticular direction/total number of genes partitioned in the block)
#partition_direction_dict = {}
for aln in Block_dict.keys():
    #aln = "Alignment 1261"
    tandem_df = pd.DataFrame(columns=["left","right"])
    #check if genes in aln_left have tandem dupicates in the tandem file
    for gene_left , gene_right in zip(Block_dict[aln]["Left"],Block_dict[aln]["Right"]):
        #check if gene left or right in tandem_df
        if gene_left in tadem_df[0].values:
            tandem_left = tadem_df[tadem_df[0] == gene_left][1].values[0]
        elif gene_left in tadem_df[1].values:
            tandem_left = tadem_df[tadem_df[1] == gene_left][0].values[0]
        else:
            tandem_left = "no"
    
        if gene_right in tadem_df[0].values:
            tandem_right = tadem_df[tadem_df[0] == gene_right][1].values[0]
        elif gene_right in tadem_df[1].values:
            tandem_right = tadem_df[tadem_df[1] == gene_right][0].values[0]
        else:
            tandem_right = "no"
        #tandem_df = tandem_df.append({"left":gene_left, "right":gene_right}, ignore_index=True)
        tandem_df = concat([tandem_df, DataFrame({"left":gene_left, "right":gene_right}, index=[0])])
        #tandem_df = tandem_df.append({"left":tandem_left, "right":tandem_right}, ignore_index=True)
        tandem_df = concat([tandem_df, DataFrame({"left":tandem_left, "right":tandem_right}, index=[0])])
    #remove rows if either left and right are "no"
    tandem_df = tandem_df[(tandem_df["left"] != "no") & (tandem_df["right"] != "no")]
    #caculate partition score
    total_gene_with_exp = 0
    total_gene_level_partition_score = 0
    for rows in tandem_df.iterrows():
        #rows = tandem_df.sample(1)
        genes = rows[1].to_list()
        
        #genes = tandem_df.iloc[3].to_list()
        
        sub_tpm_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in genes]]
        #sub_tpm_df["Dry.seeds"]
        #add very small number to every value in the df to avoid division by zero
        sub_tpm_df_mod = sub_tpm_df + 0.0000001
        #calculate fold change
        sub_tpm_df_FC = DataFrame(np.log2(sub_tpm_df_mod.iloc[0]/sub_tpm_df_mod.iloc[1]))
        #add a new column to the df with the sum of the row
        sub_tpm_df_FC["expression"] = 0
        
        #check if one of the rows is >= 2 for each column in the df sub_tpm_df
        sub_tpm_df_FC["expression"] = sub_tpm_df.apply(lambda x: 1 if x[0] >= 2 or x[1] >= 2 else 0)
        
        #add partiton score
        sub_tpm_df_FC["partition_score"] = 0
        #if 0 column > 2 and expression is 1 add 1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: 1 if x[0] >= diff_thre and x["expression"] == 1 else 0, axis=1)
        #if 1 column > 2 and expression is 1 add -1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: -1 if x[0] <= -(diff_thre) and x["expression"] == 1 else x["partition_score"], axis=1)
        
        #get total number of tissues showing the partition for syntenic gene pair
        tot_tissue_part= sub_tpm_df_FC["partition_score"].sum()
        
        # if  expression column contain 1s at 1 to total_gene_with_exp
        if sub_tpm_df_FC["expression"].sum() > 0:
            total_gene_with_exp += 1
        
        # flattening the total number of tissues showing the partition for syntenic gene pair
        if tot_tissue_part > 0:
            tot_tissue_part = 1
        elif tot_tissue_part < 0:
            tot_tissue_part = -1
        
        #add the partition score to the total_gene_level_partition_score
        total_gene_level_partition_score += tot_tissue_part
    #calculating block level partition score
    block_partition_score = total_gene_level_partition_score/total_gene_with_exp
    #add block_partition_score to partition_score_dict
    partition_score_dict[aln] = block_partition_score
#____________________________________________________________________________________#
#look at the ditribution of partition scores

#get median of partition scores for each alignment and add that in a dataframe
df_median_part = pd.DataFrame(columns=["alignmnt","block_partition_score"])
#convert partition_score_dict to a df
for aln in partition_score_dict.keys():
    df_median_part = concat([df_median_part,DataFrame({"alignmnt":aln, "block_partition_score":partition_score_dict[aln]}, index=[0])])

#show NA rows
df_median_part[df_median_part.isna().any(axis=1)]

#drop NA rows
#Na rows are dropped because the whole block does not contain any gene with expression >= 2 in any tissue
df_median_part = df_median_part.dropna()

#add a new column to df_median_part with the direction of the partition
df_median_part["block_partition_score_no_dir"] = abs(df_median_part["block_partition_score"])


#draw both plots in one figure
fig, ax = plt.subplots(1,2, figsize=(10,5))
sns.histplot(df_median_part["block_partition_score"], bins=45,element="step", fill=True, color="blue",alpha=0.5, ax=ax[0])
sns.histplot(df_median_part["block_partition_score_no_dir"], bins=45,element="step", fill=True, color="blue",alpha=0.5, ax=ax[1])
xlabel = ax[0].set_xlabel("Block level partition")
xlabel = ax[1].set_xlabel("Block level partition without direction")
#add y limit
ax[0].set_ylim(0, 55)
ax[1].set_ylim(0, 55)
#save to pdf
plt.savefig(f"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/BLOCK_level_partition_score_with_without_direction_df_{diff_thre}.pdf")
plt.close()
#remove figures from memory
plt.cla()

####################################################################################################################################################################################
#########################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
#_____________________________________________________________________________________________________________________________________#
#simulating the background distribution of partition scores
#___________________________________________________________#
#FUNCTIONS
def get_block_genes_with_tandem_dup(genome_tandems_df, Dict_blocks, alignment_nm):
    '''
    This is a function to get genes in a syntenic block with the tandem duplicated genes
    input: 
        genome_tandems_df: dataframe with tandem duplicated genes for whole genome
        Dict_blocks: dictionary with syntenic blocks information
        alignment_nm: alignment name/block to look at
    '''
    tandem_df = pd.DataFrame(columns=["left","right"])
    #check if genes in aln_left have tandem dupicates in the tandem file
    for gene_left , gene_right in zip(Dict_blocks[alignment_nm]["Left"],Dict_blocks[alignment_nm]["Right"]):
        #check if gene left or right in tandem_df
        if gene_left in genome_tandems_df[0].values:
            tandem_left = genome_tandems_df[genome_tandems_df[0] == gene_left][1].values[0]
        elif gene_left in genome_tandems_df[1].values:
            tandem_left = genome_tandems_df[genome_tandems_df[1] == gene_left][0].values[0]
        else:
            tandem_left = "no"
    
        if gene_right in genome_tandems_df[0].values:
            tandem_right = genome_tandems_df[genome_tandems_df[0] == gene_right][1].values[0]
        elif gene_right in genome_tandems_df[1].values:
            tandem_right = genome_tandems_df[genome_tandems_df[1] == gene_right][0].values[0]
        else:
            tandem_right = "no"
        #tandem_df = tandem_df.append({"left":gene_left, "right":gene_right}, ignore_index=True)
        tandem_df = concat([tandem_df, DataFrame({"left":gene_left, "right":gene_right}, index=[0])])
        #tandem_df = tandem_df.append({"left":tandem_left, "right":tandem_right}, ignore_index=True)
        tandem_df = concat([tandem_df, DataFrame({"left":tandem_left, "right":tandem_right}, index=[0])])
    #remove rows if either left and right are "no"
    tandem_df = tandem_df[(tandem_df["left"] != "no") & (tandem_df["right"] != "no")]  
    return tandem_df

#functiona to get simulated partition score
#1. Get the simulated partition score by shuffling the expression values of the tissues in syntenic gene pairs
def get_tissue_shuffel_simulated_partition_score(df_genes, df_TPM,thresh):
    '''
    This is a function to get background distribution of partition scores by shuffling expression values of the tissues in
    syntenic gene pairs
    input: df_genes: dataframe with two columns with gene pairs of syntenic block
           df_TPM: dataframe with TPM values of all the genes in the genome
           thresh : threshold to call two genes are differentially expressed based on log2(TPM_gene1/TPM_gene2) values
    '''
    total_gene_with_exp = 0
    total_gene_level_partition_score = 0
    for rows in df_genes.iterrows():
        #rows = tandem_df.sample(1)
        genes = rows[1].to_list()
        
        #genes = tandem_df.iloc[3].to_list()
        sub_tpm_df = df_TPM.loc[[re.sub(r'\.\d+$', '', x) for x in genes]]
        #draw set of random set of numbers from 0 to number of columns in the df
        rand_col = np.random.choice(range(sub_tpm_df.shape[1]),np.random.choice(range(sub_tpm_df.shape[1]),1),replace=False)
        #for rand_col swap the values in the columns
        for col in rand_col:
            sub_tpm_df.iloc[:,col] = sub_tpm_df.iloc[:,col].iloc[::-1].values
        sub_tpm_df
        #sub_tpm_df["Dry.seeds"]
        #add very small number to every value in the df to avoid division by zero
        sub_tpm_df_mod = sub_tpm_df + 0.0000001
        #calculate fold change
        sub_tpm_df_FC = DataFrame(np.log2(sub_tpm_df_mod.iloc[0]/sub_tpm_df_mod.iloc[1]))
        #add a new column to the df with the sum of the row
        sub_tpm_df_FC["expression"] = 0
        
        #check if one of the rows is >= 2 for each column in the df sub_tpm_df
        sub_tpm_df_FC["expression"] = sub_tpm_df.apply(lambda x: 1 if x[0] >= 2 or x[1] >= 2 else 0)
        
        #add partiton score
        sub_tpm_df_FC["partition_score"] = 0
        #if 0 column > 2 and expression is 1 add 1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: 1 if x[0] >= thresh and x["expression"] == 1 else 0, axis=1)
        #if 1 column > 2 and expression is 1 add -1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: -1 if x[0] <= -(thresh) and x["expression"] == 1 else x["partition_score"], axis=1)
        
        #get total number of tissues showing the partition for syntenic gene pair
        tot_tissue_part= sub_tpm_df_FC["partition_score"].sum()
        
        # if  expression column contain 1s at 1 to total_gene_with_exp
        if sub_tpm_df_FC["expression"].sum() > 0:
            total_gene_with_exp += 1
        
        # flattening the total number of tissues showing the partition for syntenic gene pair
        if tot_tissue_part > 0:
            tot_tissue_part = 1
        elif tot_tissue_part < 0:
            tot_tissue_part = -1
        
        #add the partition score to the total_gene_level_partition_score
        total_gene_level_partition_score += tot_tissue_part
    #calculating block level partition score
    return total_gene_level_partition_score/total_gene_with_exp

#2. Get the simulated partition score by shuffling the genes between the blocks
def get_block_gene_shuffel_simulated_partition_score(df_genes, df_TPM,thresh):
    '''
    This is a function to get background distribution of partition scores by shuffling genes between the blocks
    input: df_genes: dataframe with two columns with gene pairs of syntenic block
           df_TPM: dataframe with TPM values of all the genes in the genome
           trhesh : threshold to call two genes are differentially expressed based on log2(TPM_gene1/TPM_gene2) values
    '''
    total_gene_with_exp = 0
    total_gene_level_partition_score = 0
    for rows in df_genes.iterrows():
        #rows = tandem_df.sample(1)
        genes = rows[1].to_list()
        
        #genes = tandem_df.iloc[3].to_list()
        sub_tpm_df = df_TPM.loc[[re.sub(r'\.\d+$', '', x) for x in genes]]
        #pd.set_option('display.max_columns', 23)
        #sub_tpm_df
        #switch =  np.random.choice(["yes","no"],1)[0]
        if switch == "yes":
            #make row 2 row 1
            sub_tpm_df = sub_tpm_df.iloc[::-1]
        #sub_tpm_df["Dry.seeds"]
        #add very small number to every value in the df to avoid division by zero
        sub_tpm_df_mod = sub_tpm_df + 0.0000001
        #calculate fold change
        sub_tpm_df_FC = DataFrame(np.log2(sub_tpm_df_mod.iloc[0]/sub_tpm_df_mod.iloc[1]))
        #add a new column to the df with the sum of the row
        sub_tpm_df_FC["expression"] = 0
        
        #check if one of the rows is >= 2 for each column in the df sub_tpm_df
        sub_tpm_df_FC["expression"] = sub_tpm_df.apply(lambda x: 1 if x[0] >= 2 or x[1] >= 2 else 0)
        
        #add partiton score
        sub_tpm_df_FC["partition_score"] = 0
        #if 0 column > 2 and expression is 1 add 1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: 1 if x[0] >= thresh and x["expression"] == 1 else 0, axis=1)
        #if 1 column > 2 and expression is 1 add -1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: -1 if x[0] <= -(thresh) and x["expression"] == 1 else x["partition_score"], axis=1)
        
        #get total number of tissues showing the partition for syntenic gene pair
        tot_tissue_part= sub_tpm_df_FC["partition_score"].sum()
        
        # if  expression column contain 1s at 1 to total_gene_with_exp
        if sub_tpm_df_FC["expression"].sum() > 0:
            total_gene_with_exp += 1
        
        # flattening the total number of tissues showing the partition for syntenic gene pair
        if tot_tissue_part > 0:
            tot_tissue_part = 1
        elif tot_tissue_part < 0:
            tot_tissue_part = -1
        
        #add the partition score to the total_gene_level_partition_score
        total_gene_level_partition_score += tot_tissue_part
    #calculating block level partition score
    return total_gene_level_partition_score/total_gene_with_exp
#3. swapping random blocks
def swap_blocks(genome_tandems_df, Dict_blocks, alignment_nm):
    '''
    This is a function to get genes in a syntenic block with the tandem duplicated genes for left block and get a random
    block from the block pool an trim it to the size of the left block
    input: 
        genome_tandems_df: dataframe with tandem duplicated genes for whole genome
        Dict_blocks: dictionary with syntenic blocks information
        alignment_nm: alignment name/block to look at
    '''
    original_block = get_block_genes_with_tandem_dup(genome_tandems_df, Dict_blocks, alignment_nm)
    
    #pick a random key from the Dict_blocks
    random_key = np.random.choice(list(Dict_blocks.keys()),1)[0]
    #pick a random side of the block
    random_side = np.random.choice(["Left","Right"],1)[0]
    #pick a set of random genes untill len(random_genes) >= len(original_block)
    while len(Dict_blocks[random_key][random_side]) < len(original_block):
        random_key = np.random.choice(list(Dict_blocks.keys()),1)[0]
        random_side = np.random.choice(["Left","Right"],1)[0]
    genes_random = Dict_blocks[random_key][random_side]
    #pick random set of genes from the random block that equl to the number of genes in the original block
    genes_random = np.random.choice(genes_random,len(original_block),replace=False)
    if random_side == "Left":
        #replce the left block with the random block in original_block
        original_block["left"] = genes_random
    else:
        #replce the right block with the random block in original_block
        original_block["right"] = genes_random
    return original_block
#############################################################################################################################################################################
#import time
#simulation
n_iter = 100
#trshould to call two gens are differentially expressed based on log2(TPM_gene1/TPM_gene2) values
diff_thre = 2
#dict to save partition scores
sim_partition_score = []
#dict to save partition direction (this is numbers of genes partitioned in perticular direction/total number of genes partitioned in the block)
#partition_direction_dict = {}
for iter in range(n_iter):
    for aln in Block_dict.keys():
        #tandem_df = get_block_genes_with_tandem_dup(tadem_df, Block_dict, aln)
        tadem_df = swap_blocks(tadem_df, Block_dict, aln)
        #caculate partition score
        #_________________________#
        #add block_partition_score to partition_score_dict
        sim_partition_score.append(get_block_gene_shuffel_simulated_partition_score(tandem_df,TPM_df,diff_thre))
 
#remove NA values from sim_partition_score
sim_partition_score_mod = [x for x in sim_partition_score if str(x) != 'nan']

#draw the background distribution of partition scores
sns.histplot(sim_partition_score_mod, bins=45,element="step", fill=True, color="pink")
plt.xlabel("Simulated block level partition")
#plt.show()
#save to pdf
plt.savefig(f"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Simulated_bu_swapping_random_tissues_in_pairs_of_genes_BLOCK_level_partition_score_{diff_thre}.pdf")
plt.close()

####################################################################################################################################################################################
def get_actual_part_score(df_genes,df_TPM,thresh):
    '''
    This is a function to get actual partition score
    input: df_genes: dataframe with two columns with gene pairs of syntenic block
           df_TPM: dataframe with TPM values of all the genes in the genome
           thresh : threshold to call two genes are differentially expressed based on log2(TPM_gene1/TPM_gene2) values
    '''
    total_gene_with_exp = 0
    total_gene_level_partition_score = 0
    for rows in df_genes.iterrows():
        #rows = tandem_df.sample(1)
        genes = rows[1].to_list()
        
        #genes = tandem_df.iloc[3].to_list()
        
        sub_tpm_df = df_TPM.loc[[re.sub(r'\.\d+$', '', x) for x in genes]]
        #sub_tpm_df["Dry.seeds"]
        #add very small number to every value in the df to avoid division by zero
        sub_tpm_df_mod = sub_tpm_df + 0.0000001
        #calculate fold change
        sub_tpm_df_FC = DataFrame(np.log2(sub_tpm_df_mod.iloc[0]/sub_tpm_df_mod.iloc[1]))
        #add a new column to the df with the sum of the row
        sub_tpm_df_FC["expression"] = 0
        
        #check if one of the rows is >= 2 for each column in the df sub_tpm_df
        sub_tpm_df_FC["expression"] = sub_tpm_df.apply(lambda x: 1 if x[0] >= 2 or x[1] >= 2 else 0)
        
        #add partiton score
        sub_tpm_df_FC["partition_score"] = 0
        #if 0 column > 2 and expression is 1 add 1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: 1 if x[0] >= thresh and x["expression"] == 1 else 0, axis=1)
        #if 1 column > 2 and expression is 1 add -1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: -1 if x[0] <= -(thresh) and x["expression"] == 1 else x["partition_score"], axis=1)
for rows in tandem_df.iterrows():
        #rows = tandem_df.sample(1)
        genes = rows[1].to_list()
        
        #genes = tandem_df.iloc[3].to_list()
        
        sub_tpm_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in genes]]
        #sub_tpm_df["Dry.seeds"]
        #add very small number to every value in the df to avoid division by zero
        sub_tpm_df_mod = sub_tpm_df + 0.0000001
        #calculate fold change
        sub_tpm_df_FC = DataFrame(np.log2(sub_tpm_df_mod.iloc[0]/sub_tpm_df_mod.iloc[1]))
        #add a new column to the df with the sum of the row
        sub_tpm_df_FC["expression"] = 0
        
        #check if one of the rows is >= 2 for each column in the df sub_tpm_df
        sub_tpm_df_FC["expression"] = sub_tpm_df.apply(lambda x: 1 if x[0] >= 2 or x[1] >= 2 else 0)
        
        #add partiton score
        sub_tpm_df_FC["partition_score"] = 0
        #if 0 column > 2 and expression is 1 add 1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: 1 if x[0] >= diff_thre and x["expression"] == 1 else 0, axis=1)
        #if 1 column > 2 and expression is 1 add -1
        sub_tpm_df_FC["partition_score"] = sub_tpm_df_FC.apply(lambda x: -1 if x[0] <= -(diff_thre) and x["expression"] == 1 else x["partition_score"], axis=1)
        
        #get total number of tissues showing the partition for syntenic gene pair
        tot_tissue_part= sub_tpm_df_FC["partition_score"].sum()
        
        # if  expression column contain 1s at 1 to total_gene_with_exp
        if sub_tpm_df_FC["expression"].sum() > 0:
            total_gene_with_exp += 1
        
        # flattening the total number of tissues showing the partition for syntenic gene pair
        if tot_tissue_part > 0:
            tot_tissue_part = 1
        elif tot_tissue_part < 0:
            tot_tissue_part = -1
        
        #add the partition score to the total_gene_level_partition_score
        total_gene_level_partition_score += tot_tissue_part
    #calculating block level partition score
block_partition_score = total_gene_level_partition_score/total_gene_with_exp