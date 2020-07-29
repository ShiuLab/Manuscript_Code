from PIL import Image, ImageEnhance
from PIL import ImageFilter
import os,sys
if len(sys.argv) != 3:
	print('''Error: we need 3 parameters, 
	1: type of change; 
	2. degree of change; 
	3. image path.''')
type = sys.argv[1] #change type: Bright:Br, Blur:Bl, Resolution:Re, Contrast: Co 
para = sys.argv[2] #degree of change
path = sys.argv[3] #image path



fileList = os.listdir(path)
if type == "Br":
	os.mkdir(type+"_"+para)
	for f in fileList:	
		im = Image.open(path+"/"+f)
		enhancer = ImageEnhance.Brightness(im)
		enhanced_im = enhancer.enhance(float(para))
		enhanced_im.save(type+"_"+para+"/Change_Bright_%s_%s"%(para,f))
elif type == "Bl":
	os.mkdir(type+"_"+para)
	for f in fileList:	
		im = Image.open(path+"/"+f)
		blurred = im.filter(ImageFilter.GaussianBlur(radius=float(para)))
		blurred.save(type+"_"+para+"/Change_Blur_%s_%s"%(para,f))
elif type == 'Re':
	os.mkdir(type+"_"+para)
	for f in fileList:
		im = Image.open(path+"/"+f)
		im = im.resize((int(im.size[0]*float(para)), int(im.size[1]*float(para))),Image.ANTIALIAS)
		im.save(type+"_"+para+"/Change_Resolution_%s_%s"%(para,f))
elif type == 'Co':
	os.mkdir(type+"_"+para)
	for f in fileList: 
		im = Image.open(path+"/"+f)
		enhancer = ImageEnhance.Contrast(im)
		enhanced_im = enhancer.enhance(float(para))
		enhanced_im.save(type+"_"+para+"/Change_Contrast_%s_%s"%(para,f))
else:
	print('''Error: Please choice the right type: 
	Br for Bright change; 
	Bl for Blur change; 
	Re for Resolution change;
	Co for Contrast change.''')


