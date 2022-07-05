"""
This is a script you shoud rund after running batch_merge_results.py to complete the batch process of Christina's TBBM similarity pipeline
Input: Non
Output: Significantly enriched TFs for given set of kmers
Author: Thilanka Ranaweera
"""
import os
for dirs in os.listdir("./"):
    if dirs.endswith("_kmers"):
        os.chdir(dirs)
        out = open("batch_ii.sh","w")
        out.write(f'#!/bin/sh --login\n#SBATCH --time=2:00:00\n#SBATCH --ntasks=10\n#SBATCH --cpus-per-task=2\n#SBATCH --mem=20G\n#SBATCH --job-name {dirs}\n#SBATCH -e run.sh.e\n#SBATCH -o run.sh.o\ncd ./\nmodule load GCC/5.4.0-2.26 TAMO Python/2.7.10\n')
        kmer_tamo = dirs+"_no_RC.txt.tamo"
        pcc_cis = kmer_tamo+"-Athaliana_TFBM_v1.01.tm.index.direct.index.tm.dm"
        pcc_dap = kmer_tamo+"-DAP_motifs.txt.tm.dm"
        #print(f'{dirs+"_no_RC.txt.tamo"}')
        out.write(f'python /mnt/home/ranawee1/02_switchgrass_CRE_evolution/TFBM/MotifDiscovery/TAMO_scripts/pcc_merge_CC.py combine_distance_matrix_2 -t {kmer_tamo} -t2 Athaliana_TFBM_v1.01.tm.index.direct.index.tm\n')
        out.write(f'python /mnt/home/ranawee1/02_switchgrass_CRE_evolution/TFBM/MotifDiscovery/TAMO_scripts/pcc_merge_CC.py combine_distance_matrix_2 -t {kmer_tamo} -t2 DAP_motifs.txt.tm\n')
        out.write(f'python /mnt/home/ranawee1/02_switchgrass_CRE_evolution/TFBM/motif_pcc_best_match.py -t {kmer_tamo} -pcc_cis {pcc_cis} -pcc_dap {pcc_dap} -save pcre\n')
        out.write(f'python /mnt/home/ranawee1/02_switchgrass_CRE_evolution/TFBM/TFBM_pcc_significance.py -top pcre_TopHits.txt')
        out.close()
        os.system("sbatch batch_ii.sh")
        os.chdir("../")



