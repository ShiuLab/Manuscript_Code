from __future__ import print_function
import os,sys
from PIL import Image

path = "/mnt/home/mengfanr/Arabidopsis_seeds_detect/faster_rcnn_seeds_part_tiles_all_1500/detect_csv"
os.chdir(path)
for root, dirs, files in os.walk(path):
	for f in files:
		if f.startswith(sys.argv[1]+'_0_0'):
			f1=f
		elif f.startswith(sys.argv[1]+'_1_0'):
            		f2=f
		elif f.startswith(sys.argv[1]+'_0_1'):
			f3=f
		elif f.startswith(sys.argv[1]+'_1_1'):
			f4=f
out=open(sys.argv[1]+"_all.csv","a")
file1 = open(f1,"r").readlines()
for line  in file1:
	out.write("%s\n"%line.strip())

file2 = open(f2,"r").readlines()
for line  in file2:
	li = line.strip().split(",")
	name = li[0]
	width = li[1]
	height = li[2]
	cls = li[3]
	#xmin,ymin,xmax,ymax
	xmin = int(li[4])
	ymin = int(li[5])+1450
	xmax = int(li[6])
	ymax = int(li[7])+1450
	
        out.write("%s,%s,%s,%s,%d,%d,%d,%d\n"%(name,width,height,cls,xmin,ymin,xmax,ymax))

file3 = open(f3,"r").readlines()
for line  in file3:
	li = line.strip().split(",")
        name = li[0]
        width = li[1]
        height = li[2]
        cls = li[3]
        #xmin,ymin,xmax,ymax
        xmin = int(li[4])+1450
        ymin = int(li[5])
        xmax = int(li[6])+1450
        ymax = int(li[7])


        out.write("%s,%s,%s,%s,%d,%d,%d,%d\n"%(name,width,height,cls,xmin,ymin,xmax,ymax))

file4 = open(f4,"r").readlines()
for line  in file4:
        li = line.strip().split(",")
        name = li[0]
        width = li[1]
        height = li[2]
        cls = li[3]
        #xmin,ymin,xmax,ymax
        xmin = int(li[4])+1450
        ymin = int(li[5])+1450
        xmax = int(li[6])+1450
        ymax = int(li[7])+1450


        out.write("%s,%s,%s,%s,%d,%d,%d,%d\n"%(name,width,height,cls,xmin,ymin,xmax,ymax))




