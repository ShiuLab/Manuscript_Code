'''
Author LT
This is a code to check if boostrap tress have finished running. For unfinished trees shell scripts and aligned fas files will be moved to incomplete_trees dir
input: number of groups. Should always be a int.
To run this tree move all the files (if you have more than 1k groups you might have sepRted them into couple folders) into singular script
'''
import os,sys

num_groups = int(sys.argv[1]) # number of syntenic homelog groups
tree_dir = sys.argv[2] # this is the directory where the trees are located

os.chdir(tree_dir)

if "incomplete_trees" not in os.listdir("./"):
    os.mkdir("incomplete_trees")

for num in range(1,(num_groups+1),1):
    #print(f"RAxML_bipartitions.Group{num}")
    tree = [i for i in os.listdir("./") if i.startswith(f"RAxML_bipartitions.Group{num}")]
    fas = [j for j in os.listdir("./") if (j.startswith(f"Group{num}") and (not j.endswith("_aligned.fas")) and (j.endswith(".fas")))]
    sh = [k for k in os.listdir("./") if (k.startswith(f"Group{num}") and (k.endswith(".sh")))]
    if len(tree) == 0:
        os.system(f'mv {fas[0]} incomplete_trees')
        os.system(f'mv {sh[0]} incomplete_trees')