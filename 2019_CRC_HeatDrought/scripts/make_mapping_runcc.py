#make_mapping_runcc.py

import os, sys

for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-fasta':    #Path to directory with files
    path = sys.argv[i+1]
  if sys.argv[i] == '-pwm':
    pwm = sys.argv[i+1]
  if sys.argv[i] == '-include':
    include = sys.argv[i+1]

out_name = "runcc_pwm_mapping.txt"
out = open(out_name, 'w')


for f in os.listdir(path):
      if f.startswith(".") or not ".fa" in f:
        pass

      else:
      	out.write("python /mnt/home/azodichr/GitHub/Utilities/Map_PWM_biopy.py -m %s -fasta %s -include %s -out %s\n"  % (pwm, path+f, include, f+"_hits.txt"))


#os.system("python ~shius/codes/qsub_hpc.py -f queue -u azodichr -c %s -w 239 -m 4 -n 900 -wd /mnt/home/azodichr/01_DualStress_At/14_LogicClusters/17_Hybrid_Features/00_CISbp_mapping" % (file1, name, name2))
