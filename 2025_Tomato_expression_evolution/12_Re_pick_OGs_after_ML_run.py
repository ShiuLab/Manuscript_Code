'''
This code pick new OGs for the syntenic homeolg blocks whre previous OGs clads inside the ingroup during
ML tree build.
Input: 1. Alignmnt guide file from previous run /This is the file you get after running 04_check_tree_topology.py/
       2. All agains all blast result file you use to get the OGs
Output:1. Modified Alignmnt guide file. This code only modify the "out_group" coulmn. Then when you run 04_check_tree_topology.py
"topology" coulmn will be updated
       2. Alignmnt guide like file which you will use to run 01_get_fasta_for_gene_lists.py to get fasta
       //Making 2nd output becase we cannot use Modified Alignmnt guide file to get the fasta file the way code is setup//
        
'''
import pandas as pd
import itertools
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import numpy as np
import copy
import argparse
import os

#functions
def blast_df_lst(list_gene,df_blast):
    '''This is a function to output sorted blast df for list of genes
    Input: list of genes to subset from dataframe based on column 0 and column 1
            all against all blastp dataframe in output format 6
    Output: Sorted sub dataframe for list of genes
    '''
    sub_df = df_blast[(df_blast[0].isin(list_gene))| (df_blast[1].isin(list_gene))]
    soted_sub_df = sub_df.sort_values(by = 10, ascending = False)
    return soted_sub_df

def get_og_lst(sorted_df_blast,list_gene, num_og):
    '''
    This is a function to get the outgrops for list of genes
    Input: Sorted sub blast dataframe for list of genes
           list of genes: genes belonged to syntanic gene block and tandem duplicates
           number of out groups you want to specify
    Output: Sigle outgroup for corresponding list of genes
    '''
    og_list = sorted_df_blast[1].iloc[:int(num_og)].astype('str').tolist()
    og_list.extend(sorted_df_blast[0].iloc[:int(num_og)].astype('str').tolist())
    og_final = [x for x in og_list if x not in list_gene]
    return list(set(og_final))

def E_val_tree(blast_df,gene_list_test, out_groups_ini):
    '''
    This is a function to construct a distance tree for a set of syntanic homeologs based on the blast Evalues
    This function first look for potential OG for set of syntanic homelogs. Then based on 1/-log10(evals) it construct a UPGMA tree
    input:
        1. blast_df: all agains all blast dataframe for species of interset
        2. KS_gene_dict: Dictinory saving all the syntanic paralogs \both homelogs and tandem duplicates\ per gene \sa the key\
        3. key_dict: gene of interst as a key in our KS dictinory
    output: UPGMA tree built 
    '''


    
    gene_list_in_pre_og = copy.deepcopy(gene_list_test)
    gene_list_in_pre_og.extend(out_groups_ini)
    #print(gene_list_test)
    #KS_gene_dict["Solyc03g026373.1.1"]
    sub_blast_df = blast_df_lst(gene_list_in_pre_og,blast_df) # get the subset of blast dataframe corresponding to all the items in the ingroup
    soted_sub_blast_df = sub_blast_df.sort_values(by = 10, ascending = False) #sorting sub df by e val
    soted_sub_blast_df.reset_index(drop=True, inplace=True) #resetting idexes of the sub df as it might have old index from initial dataframe
    

    combinations = list(itertools.product(gene_list_in_pre_og, repeat=2))
    combinations
    index_list = []
    for comb_i in combinations:
        #print(comb_i[0], comb_i[1])
        index_list.extend(soted_sub_blast_df[((soted_sub_blast_df[0]== comb_i[0])&(soted_sub_blast_df[1]== comb_i[1])|(soted_sub_blast_df[1]== comb_i[0])&(soted_sub_blast_df[0]== comb_i[1]))].index.tolist())

    #print(min(index_list))    
    ingroup_blast_df = soted_sub_blast_df[min(index_list):]    


    #ingroup_blast_df
    og_final = set()
    col_0_list = set(list(ingroup_blast_df[0]))
    col_1_list = set(list(ingroup_blast_df[1]))
    in_eval_list = list(ingroup_blast_df[10])
    #print(in_eval_list)
    in_group_list = list(col_0_list.union(col_1_list))
    if min(index_list) != 0: #if minimum index of the ingroup combinations are not
        for index_i in reversed(range(min(index_list))):
            switch_deep_blast = 0 #switch to check if we check for potential new OGs usin newly subsetted df
            #print(index_i)
            og_group = [soted_sub_blast_df[0].iloc[index_i],soted_sub_blast_df[1].iloc[index_i]]
            e_og = soted_sub_blast_df[10].iloc[index_i]
            og_i = [x for x in  og_group if x not in gene_list_in_pre_og]
            if (og_i[0] not in in_group_list) and (e_og > max(in_eval_list)):
                #print(soted_sub_blast_df.iloc[index_i])
                #now for the ingroup check if the e vales between (ingroup_genes - outgroup) combinations are higher than ingroup-ingroup evals
                in_og_evals_list = []
                for i_in_og in list(itertools.product(gene_list_in_pre_og, [og_i[0]])):
                    #e =  UGMATree.distance(i_in_out[0],  i_in_out[1]) 
                    in_og_evals = list(soted_sub_blast_df[((soted_sub_blast_df[0]== i_in_og[0])&(soted_sub_blast_df[1]== i_in_og[1])|(soted_sub_blast_df[1]== i_in_og[0])&(soted_sub_blast_df[0]== i_in_og[1]))][10])
                    #print(og_i[0],in_og_evals)
                    in_og_evals_list.extend(in_og_evals)
                if min(in_og_evals_list) >  max(in_eval_list):
                    og_final.add(og_i[0])
                #print(index_i)
            if len(og_final)  == 2:
                #in_group = ",".join(gene_list)
                #og = ",".join(list(og_final))
                #out_aln.write(f'{in_group}\t{og}\n')
                #print("out_of_loop")
                print(f'ingroups {gene_list_test} outgroups are {og_final}')
                break
        if len(og_final) < 2: #if we did not find two out groups
            #if we only found a outgroup we have to subset a new dataframe of blast with corresponding hits with the OG
            #Then check for potential new OGs usin newly subsetted df
            if len(og_final) == 1: 
                switch_deep_blast = 1
                new_lst = copy.deepcopy(gene_list_in_pre_og)
                new_lst.extend(list(og_final)) #adding new list of genes with found OG
                extned_blast_df_with_og = blast_df_lst(new_lst,blast_df) #new dataframe of blast corresponding to the all the ingroups + found OG
                soted_extned_blast_df_with_og = extned_blast_df_with_og.sort_values(by = 10, ascending = False) #sorting sub df by e val
                soted_extned_blast_df_with_og.reset_index(drop=True, inplace=True)
                index_list_og = [] #new index list for dataframe with OG
                for comb_i in combinations: # we are still using the same combinations of ingroups. IMPoRTANT that we are not using OG combinations here
                    index_list_og.extend(soted_extned_blast_df_with_og[((soted_extned_blast_df_with_og[0]== comb_i[0])&(soted_extned_blast_df_with_og[1]== comb_i[1])|(soted_extned_blast_df_with_og[1]== comb_i[0])&(soted_extned_blast_df_with_og[0]== comb_i[1]))].index.tolist()) #check both column 1 and column 2
                for index_j in reversed(range(min(index_list_og))):
                    #print(index_i)
                    og_group = [soted_extned_blast_df_with_og[0].iloc[index_j],soted_extned_blast_df_with_og[1].iloc[index_j]] #new oG group with new Blast df
                    og_j = [x for x in  og_group if x not in gene_list_in_pre_og] #potential OG that are not in ingroup gene list
                    if og_j[0] not in in_group_list:
                        #now for the ingroup check if the e vales between (ingroup_genes - outgroup) combinations are higher than ingroup-ingroup evals
                        in_og_evals_list = []
                        for i_in_og in list(itertools.product(new_lst, [og_j[0]])):
                            #e =  UGMATree.distance(i_in_out[0],  i_in_out[1]) 
                            in_og_evals = list(soted_extned_blast_df_with_og[((soted_extned_blast_df_with_og[0]== i_in_og[0])&(soted_extned_blast_df_with_og[1]== i_in_og[1])|(soted_extned_blast_df_with_og[1]== i_in_og[0])&(soted_extned_blast_df_with_og[0]== i_in_og[1]))][10])
                            #print(og_i[0],in_og_evals)
                            in_og_evals_list.extend(in_og_evals)
                        if min(in_og_evals_list) >  max(in_eval_list):
                            og_final.add(og_j[0])
                    if len(og_final)  == 2:
                        #in_group = ",".join(gene_list)
                        #og = ",".join(list(og_final))
                        #out_aln.write(f'{in_group}\t{og}\n')
                        #counter_1 += 1
                        #print("out_of_loop")
                        print(f'ingroups {gene_list_in_pre_og} outgroups are {og_final}')
                        break
            else:
                print(f'Warning 01_2#### need to find a better outgroup for {gene_list_test}; curring outgroup found{og_final}')
    else:
        print(f'Warning 01_1#### Cannot find outgroup for {gene_list_test}')

    #making the UPGMA tree
    in_out_list = copy.deepcopy(gene_list_test)
    in_out_list.extend(list(og_final))
    #print(gene_list_test)
    #Next is to calculate distance UPGMA treee
    df_dist = pd.DataFrame({'gene_1':[],
                        'gene_2':[],
                       '1/log(e)':[]})
    #get all combinations between list items
    combs_fr_mat = list(itertools.product(in_out_list, repeat=2))
    all_gene_list = copy.deepcopy(gene_list_test)
    all_gene_list.extend(list(og_final)) 
    all_gene_blast_df = blast_df_lst(all_gene_list,blast_df)
    for comb_i in combs_fr_mat:
        #print(comb_i[0])
        comb_df = all_gene_blast_df[((all_gene_blast_df[0]== comb_i[0])&(all_gene_blast_df[1]== comb_i[1])|(all_gene_blast_df[1]== comb_i[0])&(all_gene_blast_df[0]== comb_i[1]))]
        #subsetting using blast_df is time consuming. need to use already subsetted dfs
        #if switch_deep_blast == 1:
            #comb_df = extned_blast_df_with_og[((extned_blast_df_with_og[0]== comb_i[0])&(extned_blast_df_with_og[1]== comb_i[1])|(extned_blast_df_with_og[1]== comb_i[0])&(extned_blast_df_with_og[0]== comb_i[1]))]
        #else:
            #comb_df = soted_sub_blast_df[((soted_sub_blast_df[0]== comb_i[0])&(soted_sub_blast_df[1]== comb_i[1])|(soted_sub_blast_df[1]== comb_i[0])&(soted_sub_blast_df[0]== comb_i[1]))]

        #print(len(comb_df))
        if len(comb_df) >= 1: # when more than one blast hits for corresponding gene pair
            e_val = np.mean(list(comb_df[10]))
            if e_val < 0.1:
                #print(comb_i,e_val,-(np.log10(e_val)))
                df_cont = pd.DataFrame({'gene_1':[comb_i[0]],
                            'gene_2':[comb_i[1]],
                           '1/log(e)':[1/-(np.log10(e_val))]})
                df_dist = pd.concat([df_dist, df_cont]).reset_index(drop=True)
            else:
                df_cont = pd.DataFrame({'gene_1':[comb_i[0]],
                            'gene_2':[comb_i[1]],
                           '1/log(e)':[1]})
                df_dist = pd.concat([df_dist, df_cont]).reset_index(drop=True)
        elif len(comb_df) == 0: #gene pair has no blast hit
            df_cont = pd.DataFrame({'gene_1':[comb_i[0]],
                        'gene_2':[comb_i[1]],
                       '1/log(e)':[1]})
            df_dist = pd.concat([df_dist, df_cont]).reset_index(drop=True)

    df_dist_mat = df_dist.pivot(index = 'gene_1', columns = 'gene_2', values='1/log(e)').fillna(0) #converting the table to a matrix format
    df_dist_mat = df_dist_mat.where(np.tril(np.ones(df_dist_mat.shape)).astype(np.bool_)) #getting lower tranguler matrix
    #lower tranguler matrix will have NaNs for repeting values
    #next step is to remove repeting values.
    df_dist_mat_list = np.array(df_dist_mat.replace(np.nan, '', regex=True)).tolist() #replace NaN values with empty strings
    #Next to remove empty strings from list of lists of out matrix
    matrix = [[ch for ch in x  if ch != ''] for x in df_dist_mat_list] #this will be out matrix object
    names = list(df_dist_mat) #names of the columns of the matrix
    dm = Phylo.TreeConstruction.DistanceMatrix(names,matrix)
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()
    # Construct the phlyogenetic tree using UPGMA algorithm
    UGMATree = constructor.upgma(dm)
    #Phylo.draw(UGMATree)   

    counter = 0
    for i_in_out in list(itertools.product(gene_list_test, list(og_final))):
        distance_in_out =  UGMATree.distance(i_in_out[0],  i_in_out[1])
        for i_in_in in list(itertools.combinations(gene_list_test,2)):
            distance_in_in = UGMATree.distance(i_in_in[0],  i_in_in[1])
            if distance_in_out < distance_in_in:
                print(f'ERROR##############################################################{i_in_out}')
                counter += 1
    if counter == 0:
        print("tree outgroups are phrased correctly")
    if (counter == 0) and (len(og_final) == 2):
        return list(og_final)
    else:
        return "non"

############################################################################################################
#main
parser = argparse.ArgumentParser(description='This code pick new OGs for the syntenic homeolg blocks whre previous OGs clads inside the ingroup during ML tree build.')
parser.add_argument('-m', '--mod_df',required=True, help="Alignmnt guide file from previous run /This is the file you get after running 04_check_tree_topology.py/")
parser.add_argument('-b', '--blast_df', required=True, help="All agains all blast result file you use to get the OGs")
parser.add_argument('-r', '--run_number', required=True, help="Run number of the previous run")
parser.add_argument('-o', '--out_path', required=False,default="./", help="Path of the directory to save the output files. Default is the current working directory.")
args = parser.parse_args()
mod_df_nm = args.mod_df #"Run_6_modified_alignmnts_guide_with_key_info_with_groups.txt"
blast_df_nm =  args.blast_df #"../ITAG4.0_Tomota.1e-0_maxseq500.blast"
run_number = args.run_number

mod_df = pd.read_table(mod_df_nm)
blast_df = pd.read_csv(blast_df_nm,header=None,sep ="\t")
mod_df_not_abide = mod_df[mod_df.toplogy == "not_abide"]
#output
out_aln = open(os.path.join(args.out_path,f"rerun_{run_number}_alignmnts_guide_with_key_info.txt"),"w")
out_aln.write(f'in_group\tkey\tout_group\tGroup\n')
for ind in mod_df_not_abide.index:
    #print(mod_df_not_abide.loc[ind][0])
    ingroup = mod_df_not_abide.loc[ind][0].split(",")
    pre_og = mod_df_not_abide.loc[ind][1].split(",")
    group = mod_df_not_abide.loc[ind][2]
    if E_val_tree(blast_df,ingroup,pre_og) != "non":
        print(group)
        out_group = E_val_tree(blast_df,ingroup,pre_og)
        og = ",".join(out_group)
        out_aln.write(f'{mod_df_not_abide.loc[ind][0]}\t{mod_df_not_abide.loc[ind][1]}\t{og}\t{mod_df_not_abide.loc[ind][2]}\n')
        mod_df.loc[mod_df.Group == group, "out_group"] = og
out_aln.close()
out_nm = f'Run_{run_number}_modified_alignmnts_guide_with_key_info_with_groups.txt'
mod_df.to_csv(os.path.join(args.out_path,out_nm),sep="\t",index = False)
