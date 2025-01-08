#this is a script t0 calculate the coevolution of expression across tissue pairs
############################################################################################################
import pandas as pd
import os,sys
import re
import argparse
import itertools

#___________________________________________________________________________________________________________#
#arguments
parser = argparse.ArgumentParser(description='This script is to calculate the coevolution of expression across tissue pairs')
parser.add_argument('-tpm', '--tpm_file', type=str, help='TPM table for all tissues', required=True)
parser.add_argument('-group', '--group_type_file', type=str, help='Gain loss group types for each syntenic group', required=True)
parser.add_argument('-sum', '--sum_dir', type=str, help='Directory where summery of ASR for each tissue type is stored', required=True)
parser.add_argument('-type', '--gain_loss_type', type=str, choices= ["ii","v"],help='Gain loss type to consider. only checking between pairs that has either 1 loss (ii) or 1 gain(v)', required=False)
parser.add_argument('-t', '--tresh', type=int, help='number of the tissues less than atleast one gene should express to excule broadly expressed genes', required=True)
parser.add_argument('-o', '--outdir', type=str, help='directory to save output', required=True)
args = parser.parse_args()

#paths
tpm_file = args.tpm_file
#tpm_file = '/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/TPM_table_non_T_tissues.txt'
group_type_file = args.group_type_file
#group_type_file = "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Analysis_of_exp_coevolution/gain_loss_group_types_across_tissues_for_each_syentenic_group.txt"
#directory where summery of ASR for each tissue type is stored
sum_dir = args.sum_dir
#sum_dir = "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates"
gain_loss_type = args.gain_loss_type
#gain_loss_type = "ii"
#___________________________________________________________________________________________________________#
#data
group_types =  pd.read_csv(group_type_file, sep='\t',index_col=0)
df_tpm = pd.read_csv(tpm_file, sep="\t", index_col=0)
#remove 'Mature.anthers' column from df_tpm
df_tpm = df_tpm.drop(columns=['Mature.anthers'])
#get all the tissue pairs combinations without replacemnet 
tissue_pairs_all = list(itertools.combinations(df_tpm.columns, 2))

#___________________________________________________________________________________________________________#
depndent_indenpent_table = pd.DataFrame(columns=["tissue_1","tissue_2","number_dependent","number_depndent_expected","number_independent","number_independent_expected"])

for tissue_pairs in tissue_pairs_all:
    #print(list(tissue_pairs))
    tissue_pairs =  list(tissue_pairs)
    #tissue_pairs = ["ZZTrichome","Seedling.root.apex"]
    #tissue_pairs = ["Seedling.root.apex","ZZTrichome"]
    tissue_pair_group_types = group_types[tissue_pairs]

    #get the rows where both columns have the gain_loss_type
    sub_df_gain_loss = tissue_pair_group_types[(tissue_pair_group_types[tissue_pairs[0]]==gain_loss_type) & (tissue_pair_group_types[tissue_pairs[1]]==gain_loss_type)]

    #tpm_Dfs
    df_tpm_tissue = df_tpm[tissue_pairs]
    #get the syntenic groups where both tissues are dependent interms on expression
    number_dependent = 0
    #number independent groups
    number_independent = 0
    #table_to_save_depende_independent values for tissue pair

    for index in sub_df_gain_loss.iterrows():
        #print(index[0])
        #index = ["Group570"]
        tissue_1_path = os.path.join(sum_dir,tissue_pairs[0])
        tissue_2_path = os.path.join(sum_dir,tissue_pairs[1])
        #get the file starts with index in tissue_1_path
        group_file_1 = [f for f in os.listdir(tissue_1_path) if f.startswith(str(index[0])+"_")]
        #get the file starts with index in tissue_2_path
        group_file_2 = [f for f in os.listdir(tissue_2_path) if f.startswith(str(index[0])+"_")]
        if (len(group_file_1) != 0 ) & (len(group_file_2) != 0):
            #open_two files
            file_1 = pd.read_csv(os.path.join(tissue_1_path,group_file_1[0]), sep='\t')
            file_2 = pd.read_csv(os.path.join(tissue_2_path,group_file_2[0]), sep='\t')
            #remove ogs from files
            outgroups =  group_file_1[0].replace(".txt","").split("_")[1:3]
            
            #remove the rows if there is atleast on outgroup in node_name column
            file_1_in_group = file_1[~file_1['node_name'].str.contains('|'.join(outgroups))]
            file_2_in_group = file_2[~file_2['node_name'].str.contains('|'.join(outgroups))]

            #get in group genes
            in_group_genes = list(set([item for sublist in [x.split(",") for x in file_1_in_group['node_name'].tolist()] for item in sublist]))
            
            #getting TPM values for all the tissues for in_group_genes
            df_tpm_all_tissue = df_tpm.loc[[re.sub(r'\.\d+$', '', x) for x in in_group_genes]]
            #check if atlest one gene is expressed in more than 10 tissues
            if (df_tpm_all_tissue > 2).sum(axis=1).max() <= args.tresh:
            
                #get the TPM values for the in_group_genes
                df_tpm_in_group = df_tpm_tissue.loc[[re.sub(r'\.\d+$', '', x) for x in in_group_genes]]
                #add a new column to the df_tpm_in_group to store the node names and add ingroup genes if they match with the index
                df_tpm_in_group['node_name'] = [x for x in in_group_genes if [y in re.sub(r'\.\d+$', '', y) for y in df_tpm_in_group.index]]
                df_tpm_in_group.reset_index(drop=True, inplace=True)
                
                #add node names and respective tissue TPM values to file_1_in_group and file_2_in_group
                df_ingroup_tissue_1 =  df_tpm_in_group[["node_name",tissue_pairs[0]]]
                df_ingroup_tissue_1.rename(columns={tissue_pairs[0]: 'ancestral_val'}, inplace=True)
                #concatinate file_1_in_group and df_ingroup_tissue_1 row wise
                file_1_in_group = pd.concat([file_1_in_group,df_ingroup_tissue_1], axis=0)
                
                
                df_ingroup_tissue_2 =  df_tpm_in_group[["node_name",tissue_pairs[1]]]
                df_ingroup_tissue_2.rename(columns={tissue_pairs[1]: 'ancestral_val'}, inplace=True)
                #concatinate file_2_in_group and df_ingroup_tissue_2 row wise
                file_2_in_group = pd.concat([file_2_in_group,df_ingroup_tissue_2], axis=0)

                #make the table binary
                file_1_in_group['ancestral_val'] = file_1_in_group['ancestral_val'].apply(lambda x: 0 if x < 2 else 1)
                file_2_in_group['ancestral_val'] = file_2_in_group['ancestral_val'].apply(lambda x: 0 if x < 2 else 1)
                
                #check if both tissues are dependent
                #check for each row ancestral_val column are equal in both files
                if file_1_in_group['ancestral_val'].equals(file_2_in_group['ancestral_val']):
                    number_dependent += 1
                else:
                    number_independent += 1
                    print(index[0])
        number_expected = round((number_dependent + number_independent)/2)
        #fill the table using concat
    depndent_indenpent_table = pd.concat([depndent_indenpent_table,pd.DataFrame([[tissue_pairs[0],tissue_pairs[1],number_dependent,number_expected,number_independent,number_expected]], columns=["tissue_1","tissue_2","number_dependent","number_depndent_expected","number_independent","number_independent_expected"])], axis=0)

    #save the table
    depndent_indenpent_table.to_csv(os.path.join(args.outdir,f"dependent_independent_table_across_tissue_pairs_for_{args.gain_loss_type}_excluding_tissues_more_than_{args.tresh}.txt"), sep='\t', index=False)
        
        
    
