import os, sys
import time
from PIL import Image,ImageDraw
##absolute path of scanned images 
path = "/*/*/*"
os.chdir(path)
all_images = [f for f in os.listdir(path ) if f.endswith("jpg")]

print(all_images)
##create directory for splited images
if not os.path.exists(path+"\\single_plate"):
	os.mkdir("single_plate")
else:
        #if single_plate directory exist, backup and recreate a new one
        os.rename("single_plate","single_plate_"+time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime()))
        os.mkdir("single_plate")
Image.MAX_IMAGE_PIXELS = 10000000
for f in all_images:
        im=Image.open(path+"\\"+f)
        width=im.size[0]
        height=im.size[1]
        ##height
        for i in range(4):
                ##width
                for  j in range(3):
                        start_width=710
                        start_height=470
                        region=((start_width+   (370+2900)*(j)),(start_height+(500+2900)*(i)),(start_width+(370*j)+(2900)*(j+1)),(start_height+(500*i)+(2900)*(i+1)))
                        cropImg = im.crop(region)
                        cropImg.save("single_plate\\"+f+"_"+str(i)+"_"+str(j)+".jpg")
		