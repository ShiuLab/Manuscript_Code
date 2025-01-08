import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
#_______________________________________________________________________________________________________________________#
#getting the alignmnts with median KS values between 0.5871012 and 0.89355>> evolved from the most recent GT event
ks_df = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/meadian_ks_for_syntanic_block.txt")
sub_ks_df = ks_df[ks_df["median_KS"] > 0.5871012]
sub_ks_df =  sub_ks_df[sub_ks_df["median_KS"] < 0.89355] 
algnment_list = list(sub_ks_df["block"])

##_______________________________________________________________________________________________________________________#
ks_file = "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/ks_ka_tomato_paralogs_MCScanx_default"
tadem_file = "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/ITAG4.0_Tomota.tandem"
tadem_df =  pd.read_csv(tadem_file, header=None)
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

len(Block_dict.keys())
#draw a plot with lenth of the blocks frquency
block_len_list = []
for aln_name, blocks in Block_dict.items():
    block_len_list.append(len(blocks["Left"]))
df_block_len = pd.DataFrame({"Block_len":block_len_list})
#sort the dataframe
df_block_len = df_block_len.sort_values(by="Block_len")
print(df_block_len.to_string())
#plot the histogram
plt.hist(df_block_len["Block_len"], bins=100)
plt.xlabel("Block length")
plt.ylabel("Frequency")
plt.savefig("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Block_length_distribution.pdf")
plt.close()


out = open("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Synteny_blocks_genes_alignment_intersection.txt","w")
out.write("Aln1\tAln2\tLeft_intersection\tRight_intersection\tLeft_right_intersection\tRight_left_intersection\n")
for aln_name_i, blocks_i in Block_dict.items():
    for aln_name_j, blocks_j in Block_dict.items():
        if aln_name_j != aln_name_i:
            #get intersection with all combinations of left and right gene sets
            #print(aln_name_i, aln_name_j)
            left_i = set(blocks_i["Left"])
            left_j = set(blocks_j["Left"])
            right_i = set(blocks_i["Right"])
            right_j = set(blocks_j["Right"])
            left_intersection = left_i.intersection(left_j)
            right_intersection = right_i.intersection(right_j)
            left_right_intersection = left_i.intersection(right_j)
            right_left_intersection = right_i.intersection(left_j)
            #print(aln_name_i, aln_name_j, len(left_intersection), len(right_intersection), len(left_right_intersection), len(right_left_intersection))
            out.write(aln_name_i + "\t" + aln_name_j + "\t" + str((len(left_intersection)/len(left_i))*100) + "\t" + str((len(right_intersection)/len(right_i))*100) + "\t" + str((len(left_right_intersection)/len(left_i))*100) + "\t" + str((len(right_left_intersection)/len(right_i))*100) + "\n") 
out.close()

#get the alighnmnts that does not have any gene intersection with any other alignment
aln_intersect_df = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Synteny_blocks_genes_alignment_intersection.txt")
#change the threshold from 0 to 100 and get len(aln_with_two_blocks) in a dataframe
threshold_list = []
aln_count_list = []
for threshold in range(0,100,1):
    unique_aln_1 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] <= threshold) & (aln_intersect_df["Right_intersection"] <= threshold) & (aln_intersect_df["Left_right_intersection"] <= threshold) & (aln_intersect_df["Right_left_intersection"] <= threshold)]["Aln1"].unique())
    unique_aln_2 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] <= threshold) & (aln_intersect_df["Right_intersection"] <= threshold) & (aln_intersect_df["Left_right_intersection"] <= threshold) & (aln_intersect_df["Right_left_intersection"] <= threshold)]["Aln2"].unique())
    aln_with_two_blocks = list(set(unique_aln_1 + unique_aln_2))
    mutiple_aln_1 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] > threshold) | (aln_intersect_df["Right_intersection"] > threshold) | (aln_intersect_df["Left_right_intersection"] > threshold) | (aln_intersect_df["Right_left_intersection"] > threshold)]["Aln1"].unique())
    muitiple_aln_2 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] > threshold) | (aln_intersect_df["Right_intersection"] > threshold) | (aln_intersect_df["Left_right_intersection"] > threshold) | (aln_intersect_df["Right_left_intersection"] > threshold)]["Aln2"].unique())
    multiple_aln_block = list(set(mutiple_aln_1 + muitiple_aln_2))
    aln_with_two_blocks_mod = list(set(aln_with_two_blocks) - set(aln_with_two_blocks).intersection(set(multiple_aln_block)))
    threshold_list.append(threshold)
    aln_count_list.append(len(aln_with_two_blocks_mod))
df_threshold = pd.DataFrame({"Threshold":threshold_list, "Aln_count":aln_count_list})
print(df_threshold.to_string())
#prepare a plot
plt.plot(df_threshold["Threshold"], df_threshold["Aln_count"])
plt.xlabel("Gene overlap percentage threshold to call a block of genes are in synteny with syntenic region with two bocks")
plt.ylabel("Number of syntenic regions with two blocks")
#plt.show()
#save as a pdf
plt.savefig("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Threshold_vs_aln_count.pdf")
plt.close()
#_______________________________________________________________________________________________________________________#
final_trenchold = 40
unique_aln_1 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] <= final_trenchold) & (aln_intersect_df["Right_intersection"] <= final_trenchold) & (aln_intersect_df["Left_right_intersection"] <= final_trenchold) & (aln_intersect_df["Right_left_intersection"] <= final_trenchold)]["Aln1"].unique())
unique_aln_2 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] <= final_trenchold) & (aln_intersect_df["Right_intersection"] <= final_trenchold) & (aln_intersect_df["Left_right_intersection"] <= final_trenchold) & (aln_intersect_df["Right_left_intersection"] <= final_trenchold)]["Aln2"].unique())
aln_with_two_blocks = list(set(unique_aln_1 + unique_aln_2))
mutiple_aln_1 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] > final_trenchold) | (aln_intersect_df["Right_intersection"] > final_trenchold) | (aln_intersect_df["Left_right_intersection"] > final_trenchold) | (aln_intersect_df["Right_left_intersection"] > final_trenchold)]["Aln1"].unique())
muitiple_aln_2 = list(aln_intersect_df[(aln_intersect_df["Left_intersection"] > final_trenchold) | (aln_intersect_df["Right_intersection"] > final_trenchold) | (aln_intersect_df["Left_right_intersection"] > final_trenchold) | (aln_intersect_df["Right_left_intersection"] > final_trenchold)]["Aln2"].unique())
multiple_aln_block = list(set(mutiple_aln_1 + muitiple_aln_2))
aln_with_two_blocks_mod = list(set(aln_with_two_blocks) - set(aln_with_two_blocks).intersection(set(multiple_aln_block)))
len(aln_with_two_blocks_mod)
#_______________________________________________________________________________________________________________________#
#getting expression correlation for the genes in the alighnments
################################################################

TPM_df = pd.read_table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/TPM_table_non_T_tissues.txt",header=0, index_col=0)
df_block_corr = pd.DataFrame(columns=["Aln_name", "Left_corr_median", "Right_corr_median"])
for aln_name, blocks in Block_dict.items():
    if aln_name in aln_with_two_blocks_mod:
        left_genes = blocks["Left"]
        right_genes = blocks["Right"]
        left_TPM_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in left_genes]]
        right_TPM_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in right_genes]]
        
        #calculate pairewise correlation for the left and right genes
        left_corr = left_TPM_df.transpose().corr(method="pearson")
        #get the medium correlation value of the upper triangle
        left_corr_median_values = left_corr.where(np.triu(np.ones(left_corr.shape), k=1).astype(np.bool)).stack().median()
        right_corr = right_TPM_df.transpose().corr(method="pearson")
        #get the medium correlation value of the upper triangle
        right_corr_median_values = right_corr.where(np.triu(np.ones(right_corr.shape), k=1).astype(np.bool)).stack().median()
        #add largest coorrlation value to left and smallest to right
        #get largest correlation value between right_corr_median_values and left_corr_median_values
        df_block_corr = df_block_corr.append({"Aln_name":aln_name, "Left_corr_median":max(left_corr_median_values, right_corr_median_values), "Right_corr_median":min(left_corr_median_values, right_corr_median_values)}, ignore_index=True)

#draw scatter plot using seaborn
import seaborn as sns
sns.set(style="whitegrid")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(10, 10))
# Plot the total crashes
sns.set_color_codes("pastel")
sns.scatterplot(x="Right_corr_median", y="Left_corr_median", data=df_block_corr, s=50, color="red")
#make x and y axis go though 0,0
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
#make x limits and y limits equal
plt.xlim(-1, 1)
plt.ylim(-0.25, 1)


plt.ylabel("Left block median correlation")
plt.xlabel("Right block median correlation")
plt.savefig("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Left_right_block_correlation.pdf")
plt.close()


#draw Show filled contours of bivariate density plots
sns.set_style("white")
sns.set_style("ticks")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(10, 10))
# Plot the total crashes
sns.set_color_codes("pastel")
sns.kdeplot(df_block_corr["Right_corr_median"],df_block_corr["Left_corr_median"] , cmap="YlOrBr", shade=True, shade_lowest=True, cbar=True, n_levels=30)
sns.despine()

#sns.scatterplot(x="Right_corr_median", y="Left_corr_median", data=df_block_corr, s=10, color="black")
#make a legend
plt.legend(loc='upper right')
#make x and y axis go though 0,0    
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
#make x limits and y limits equal
plt.xlim(-0.6, 0.6)
plt.ylim(-1, 1)
plt.ylabel("Left block median correlation")
plt.xlabel("Right block median correlation")
plt.savefig("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Block_correlation_density.pdf")
plt.close()

#_______________________________________________________________________________________________________________________#
#getting a expected correlation value for the genes in the alighnments
#here keeping one block constant and next block randomlly selected from the alighnments
#doing a simulation on n iterations
n = 1000
df_block_corr_exp = pd.DataFrame(columns=["Aln_name", "Left_corr_median", "Right_corr_median"])
for i in range(n):
    print(i)
    for aln_name, blocks in Block_dict.items():
        if aln_name in aln_with_two_blocks_mod:
            #pick a random alignmnt for right genes
            aln_name_rand = np.random.choice(list(Block_dict.keys()))
            left_genes = blocks["Left"]
            right_genes = Block_dict[aln_name_rand]["Right"]
            left_TPM_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in left_genes]]
            right_TPM_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in right_genes]]
            #shuffle the rows of the TPM dataframe
            left_TPM_df = left_TPM_df.sample(frac=1).reset_index(drop=True)
            right_TPM_df = right_TPM_df.sample(frac=1).reset_index(drop=True)
            #calculate pairewise correlation for the left and right genes
            left_corr = left_TPM_df.transpose().corr(method="pearson")
            #get the medium correlation value of the upper triangle
            left_corr_median_values = left_corr.where(np.triu(np.ones(left_corr.shape), k=1).astype(np.bool)).stack().median()
            right_corr = right_TPM_df.transpose().corr(method="pearson")
            #get the medium correlation value of the upper triangle
            right_corr_median_values = right_corr.where(np.triu(np.ones(right_corr.shape), k=1).astype(np.bool)).stack().median()
            #add largest coorrlation value to left and smallest to right
            #get largest correlation value between right_corr_median_values and left_corr_median_values
            df_block_corr_exp = df_block_corr_exp.append({"Aln_name":aln_name, "Left_corr_median":max(left_corr_median_values, right_corr_median_values), "Right_corr_median":min(left_corr_median_values, right_corr_median_values)}, ignore_index=True)

#draw density plot for the expected correlation values
sns.set_style("white")
sns.set_style("ticks")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(10, 10))
# Plot the total crashes
sns.set_color_codes("pastel")
sns.kdeplot(df_block_corr_exp["Right_corr_median"],df_block_corr_exp["Left_corr_median"] , cmap="YlOrBr", shade=True, shade_lowest=True, cbar=True, n_levels=30)
sns.despine()
#draw scatter plot on top of the density plot
#sns.scatterplot(x="Right_corr_median", y="Left_corr_median", data=df_block_corr, s=10, color="black")

#make a legend
plt.legend(loc='upper right')
#make x and y axis go though 0,0
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
#make x limits and y limits equal
plt.xlim(-0.6, 0.6)
plt.ylim(-1, 1)

plt.ylabel("Left block median correlation")
plt.xlabel("Right block median correlation")
plt.savefig("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Block_correlation_density_rndom_distribution.pdf")
plt.close()
plt.clf()
#_______________________________________________________________________________________________________________________#
#getting a expected correlation value for the genes in the alighnments
#here keeping one block constant and next block adding random values to the TPM values
#doing a simulation on n iterations


n = 1000
df_block_corr_exp_2 = pd.DataFrame(columns=["Aln_name", "Left_corr_median", "Right_corr_median"])
for i in range(n):
    print(i)
    for aln_name, blocks in Block_dict.items():
        if aln_name in aln_with_two_blocks_mod:
            #pick a random alignmnt for right genes
            aln_name_rand = np.random.choice(list(Block_dict.keys()))
            left_genes = blocks["Left"]
            right_genes = Block_dict[aln_name_rand]["Right"]
            left_TPM_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in left_genes]]
            right_TPM_df = TPM_df.loc[[re.sub(r'\.\d+$', '', x) for x in right_genes]]
            #add random values to the right TPM dataframe
            max_tpm = right_TPM_df.max().max()
            min_tpm = right_TPM_df.min().min()
            right_TPM_df = right_TPM_df + np.random.uniform(min_tpm, max_tpm, right_TPM_df.shape)
            #shuffle the rows of the TPM dataframe
            left_TPM_df = left_TPM_df.sample(frac=1).reset_index(drop=True)
            right_TPM_df = right_TPM_df.sample(frac=1).reset_index(drop=True)
            #calculate pairewise correlation for the left and right genes
            left_corr = left_TPM_df.transpose().corr(method="pearson")
            #get the medium correlation value of the upper triangle
            left_corr_median_values = left_corr.where(np.triu(np.ones(left_corr.shape), k=1).astype(np.bool)).stack().median()
            right_corr = right_TPM_df.transpose().corr(method="pearson")
            #get the medium correlation value of the upper triangle
            right_corr_median_values = right_corr.where(np.triu(np.ones(right_corr.shape), k=1).astype(np.bool)).stack().median()
            #add largest coorrlation value to left and smallest to right
            #get largest correlation value between right_corr_median_values and left_corr_median_values
            df_block_corr_exp_2 = df_block_corr_exp_2.append({"Aln_name":aln_name, "Left_corr_median":max(left_corr_median_values, right_corr_median_values), "Right_corr_median":min(left_corr_median_values, right_corr_median_values)}, ignore_index=True)

#draw density plot for the expected correlation values
sns.set_style("white")
sns.set_style("ticks")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(10, 10))
# Plot the total crashes
sns.set_color_codes("pastel")
sns.kdeplot(df_block_corr_exp_2["Right_corr_median"],df_block_corr_exp_2["Left_corr_median"] , cmap="YlOrBr", shade=True, shade_lowest=True, cbar=True, n_levels=30)
sns.despine()
#draw scatter plot on top of the density plot
#sns.scatterplot(x="Right_corr_median", y="Left_corr_median", data=df_block_corr, s=10, color="black")

#make a legend
plt.legend(loc='upper right')
#make x and y axis go though 0,0
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
#make x limits and y limits equal
plt.xlim(-0.6, 0.6)
plt.ylim(-1, 1)
#plt.xlim(-0.2, 0.5)
#plt.ylim(-0.2, 0.8)
plt.ylabel("Left block median correlation")
plt.xlabel("Right block median correlation")
plt.savefig("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Analysis_on_syntanic_blocks/Figures/Block_correlation_density_rndom_distribution_2_shade_lowest.pdf")
plt.close()
plt.clf()