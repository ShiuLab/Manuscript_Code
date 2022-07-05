"""
    This code automates Christina's TBBM similarity pipeline for sets of pCRE lists 
    input: working directory for your main dir. This dir should include
        1.Git clone of Motif discovary pipeline
        2. scripts
            i. motif_pcc_best_match.py
            ii.merge_runcc.py
            iii. TFBM_pcc_significance.py
            iv. 3.calculate_disctance_2_file.py
        3.motif_files
            i.Dap seq index file and TAMO file
            ii. CIS_I index and TAMO file
        4. Your pCRE list named as <name>__kmers.txt
    Output: combined PWM for kmers and DPA/CIS motifs
        
"""
import os, sys

def make_input_files(wd):
    """
    This is a function takes in working dir and make dir for each kmer list and put relevent files in there
    return : nothing
    """
    for motif_list in os.listdir(wd):
        if motif_list.endswith("kmers.txt"):
            list_name = motif_list.replace(".txt","")
            #print(list_name)
            if list_name not in os.listdir("./"):
                os.mkdir(list_name)
            os.system(f'mv {motif_list} {list_name}')
            os.system(f'cp merge_runcc.py {list_name}')
            os.system(f'cp motif_files/** {list_name}')
            #os.system(f'')
            os.chdir(list_name)
            os.system(f'sed "s/\..*//" {motif_list} > {list_name + "_no_RC.txt"}')
            os.chdir("../")


workin_dir = sys.argv[1]
make_input_files(workin_dir)

for motif_dir in os.listdir("./"):
    if motif_dir.endswith("kmers"):
        os.chdir(motif_dir)
        file_no_rc = motif_dir+"_no_RC.txt"
        out = open(f'{motif_dir+".sh"}',"w")
        out.write(f'#!/bin/sh --login\n#SBATCH --time=2:00:00\n#SBATCH --ntasks=10\n#SBATCH --cpus-per-task=2\n#SBATCH --mem=20G\n#SBATCH --job-name {motif_dir}\n#SBATCH -e run.sh.e\n#SBATCH -o run.sh.o\ncd ./\nmodule load GCC/5.4.0-2.26 TAMO Python/2.7.10\n')
        out.write(f'python {workin_dir+"/MotifDiscovery/TAMO_scripts/generate_PWM.py"} {file_no_rc}\n')
        out.write(f'python {workin_dir+"/MotifDiscovery/TAMO_scripts/pcc_merge_CC.py"} create_cc_2 -t {file_no_rc+".tamo"} -t2 Athaliana_TFBM_v1.01.tm.index.direct.index.tm\n')
        out.write('mv runcc runcc_CIS.sh\n')
        out.write(f'python {workin_dir+"/MotifDiscovery/TAMO_scripts/pcc_merge_CC.py"} create_cc_2 -t {file_no_rc+".tamo"} -t2 DAP_motifs.txt.tm\n')
        out.write(f'mv runcc runcc_DAP.sh\n')
        out.write(f'python3 merge_runcc.py')
        out.close()
        os.system(f'sbatch {motif_dir+".sh"}')
        os.chdir("../")
