import sys
size = float(sys.argv[1])

inp = open("truth_ground.csv","r").readlines()
out = open("truth_ground_%s"%(sys.argv[1])+".csv","w")
for line in inp:
	if not line.startswith("file"):
		li = line.strip().split(",")
		#filename,width,height,class,xmin,ymin,xmax,ymax
		filename, width, height, cls, xmin, ymin, xmax, ymax=tuple(li)
		out.write("%s,%d,%d,%s,%d,%d,%d,%d\n"%(filename,int(int(width)*size),int(int(height)*size),cls,int(int(xmin)*size),int(int(ymin)*size),int(int(xmax)*size),int(int(ymax)*size)))
		


