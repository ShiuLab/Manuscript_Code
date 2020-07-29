from __future__ import print_function
import os,sys
from PIL import Image

path = "/mnt/home/mengfanr/Arabidopsis_seeds_detect/faster_rcnn_seeds_part_tiles_all_1500/test_image_new_half_half_results"
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
c1=f1.split(".")[-2].split("_")[1]
c2=f2.split(".")[-2].split("_")[1]
c3=f3.split(".")[-2].split("_")[1]
c4=f4.split(".")[-2].split("_")[1]
count = int(c1)+int(c2)+int(c3)+int(c4)
files = [f1,f2,f3,f4]

result = Image.new("RGB", (2900, 2900))

for index, file in enumerate(files):
  path = os.path.expanduser(file)
  img = Image.open(path)
  img.thumbnail((1450, 1450), Image.ANTIALIAS)
  x = index // 2 * 1450
  y = index % 2 * 1450
  w, h = img.size
  result.paste(img, (x, y, x + w, y + h))

result.save(os.path.expanduser(sys.argv[1]+"_merge_"+str(count)+".jpg"))
