import os, sys
import time
from PIL import Image,ImageDraw
from pathlib import Path 
##absolute path of scanned images 
path = sys.argv[1]
os.chdir(path)
all_images = [f for f in os.listdir(path ) if f.endswith("jpg")]

print(all_images)
##create directory for splited images
if not os.path.exists(path+"/split"):
	os.mkdir("split")
else:
        #if single_plate directory exist, backup and recreate a new one
        os.rename("split","split_"+time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime()))
        os.mkdir("split")
Image.MAX_IMAGE_PIXELS = 10000000000
for f in all_images:
        im=Image.open(path+"/"+f)
        width=im.size[0]
        height=im.size[1]
        ##height
        for i in range(2):
                ##width
                for  j in range(2):
                  
                        region=(j*width/2,i*height/2,(j+1)*width/2,(i+1)*height/2)
                        cropImg = im.crop(region)
                        cropImg.save(path+"/split/"+f+"_"+str(i)+"_"+str(j)+".jpg")
		
