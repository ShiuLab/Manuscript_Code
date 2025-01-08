'''
This is a code to check if the ML tree topology abides to the tree that was created using BLAST E-values.
This analysis is done as BLAST evalue criteria to pick the OGs could be too relaxed, thus OGs can be claded as
    in groups when the ML tree is formed
Input: The alignmnt guide file /with full path if its in a different directory/ created by previous code
Output: Modified alignmnt guide file
        You could also use > to print the output to a file
'''

from Bio import Phylo
import itertools
import copy, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse

#methods
def check_tree(in_group,out_group,clad_tree):
    '''
    this is a function to check if out groups are claded correctly in a given tree
    Input: 
        in_group: list of ingroups
        out_group: list of out group
        clad_tree: tree with transformed branch lengths to equal lenths as we are only intersed in tree topology
    Output: 1 or 0 depending on the OG plcement
    '''
    error_switch = 0
    for i_in_out in list(itertools.product(in_group, out_group)):
            distance_in_out =  clad_tree.distance(i_in_out[0],  i_in_out[1])
            for i_in_in in list(itertools.combinations(in_group,2)):
                distance_in_in = clad_tree.distance(i_in_in[0],  i_in_in[1])
                if distance_in_out < distance_in_in:
                    error_switch = 1
    return error_switch

#parse arguments
parser = argparse.ArgumentParser(description='Check if the ML tree topology abides to the tree that was created using BLAST E-values')
parser.add_argument('-aln', '--aln_guide',required=True, help="alignmnt guide file with full path")
parser.add_argument('-t', '--tree_dir',required=True, help="directory where the trees are located")
args = parser.parse_args()
#aln_g = sys.argv[1] # alignmnt guide file with full path

aln_guide = pd.read_table(args.aln_guide)
if 'toplogy' not in aln_guide.columns:
    aln_guide["toplogy"] = "not_abide"
    
os.chdir(args.tree_dir)
 
#counter to count trees with wrong OG plcement
count_error = 0
count_tot = 0
if not "0_correct_trees" in os.listdir("./"):
    os.mkdir("0_correct_trees")

for tree_name in os.listdir("./"):
    if tree_name.startswith("RAxML_bipartitions.Group"):
        count_tot += 1
        #tree_name =  "RAxML_bipartitions.Group1_Solyc07g062100.4.1_Solyc12g082795.1.1_RAxML_1000.out"
        tree = Phylo.read(tree_name, "newick")
        #draw tree before transforming branches
        #Phylo.draw(tree)
        ogs = tree_name.split("_")[2:4]
        group= int(tree_name.split("_")[1].replace("bipartitions.Group",""))
        tree_transfomed = copy.deepcopy(tree)
        # Modify the branch lengths
            #do this as branch lengths of original ML tree are nucleotide substitution rates. Thus some branches could be longer compared with OG
            # we only need the topology of the tree. Thus we can make all the branches to have equal lengths
        for clade in tree_transfomed.find_clades():
            clade.branch_length = 2

        #visulization of trees
        #original tree
        Phylo.draw(tree)
        # transformed brances tree
        Phylo.draw(tree_transfomed)

        #extracting tip lable info from tree object
        tip_lables = []
        lengths = []
        for terminals in tree.get_terminals():
            tip_lables.append(terminals.name)
            #lengths
            lengths.append(tree.clade.distance(terminals.name))
        ingroup = list(set(set(tip_lables) - set(ogs)))

        if check_tree(ingroup,ogs,tree_transfomed) == 1:
            print(f"###ERROR##### Out groups claded with terminal branches {ogs} \n\n\n")
            count_error += 1
        else:
            print("Correct OG plcement")
            aln_guide.loc[aln_guide.Group == group, "toplogy"] = "abide"
            os.system(f'mv {tree_name} ./0_correct_trees')
print(f"error outgroup plcement rate is {(count_error/count_tot)*100}")
#save final dataframe
aln_guide.to_csv("Modified_alignmnts_guide_with_key_info_with_groups.txt",sep="\t",index = False)


