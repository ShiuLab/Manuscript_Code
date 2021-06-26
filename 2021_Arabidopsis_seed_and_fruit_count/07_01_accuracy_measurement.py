import os,sys
import numpy as np
ground = open(sys.argv[1],"r").readlines()
detect = open(sys.argv[2],'r').readlines()
def read_boxes(boxes):
	D = {}
	for box in boxes:
		if not box.startswith("file"):
			tem = box.strip().split(',')
			image = box.strip().split(",")[0]
			if image not in D.keys():
				D[image]=[]
			D[image].append([int(tem[4]),int(tem[5]),int(tem[6]),int(tem[7])])
	return D
##xmin,ymin,xmax,ymax
def compute_iou(box1,box2):
	g_ymin, g_xmin, g_ymax, g_xmax = tuple(box1)
	d_ymin, d_xmin, d_ymax, d_xmax = tuple(box2)
	##intersection area
	xa = max(g_xmin, d_xmin)
	ya = max(g_ymin, d_ymin)
	xb = min(g_xmax, d_xmax)
	yb = min(g_ymax, d_ymax)
	intersection = max(0, xb - xa + 1) * max(0, yb - ya + 1)
	boxAArea = (g_xmax - g_xmin + 1) * (g_ymax - g_ymin + 1)
	boxBArea = (d_xmax - d_xmin + 1) * (d_ymax - d_ymin + 1)
	return intersection / float(boxAArea + boxBArea - intersection)

groundtrue = read_boxes(ground)
detected = read_boxes(detect)

iou_D = {}
out = open(sys.argv[2].split(".csv")[0]+"_accuracy.csv","w")
out.write("Name,True,Predection,TP,FP,FN,Precision,Recall,Accuracy,F1\n")
results = []
for image in detected:
	matches = []
	for i in range(len(detected[image])):
		for j in range(len(groundtrue[image])):
			iou = compute_iou(groundtrue[image][j],detected[image][i])
			if iou >= 0.5:
				matches.append([i, j, iou])
	matches = np.array(matches)
	TP = len(np.unique(matches[:,1]))
	Trues = len(groundtrue[image])
	Detected = len(detected[image])
	FP = Detected - TP
	FN = Trues - TP
	Precision = '{:.6f}'.format(float(TP/(TP+FP)))
	Recall = '{:.6f}'.format(float(TP/(TP+FN)))
	Accuracy = '{:.6f}'.format(float(TP/(TP+FN+FP)))
	F1 = '{:.6f}'.format(float(2*(float(Precision)*float(Recall)/(float(Precision)+float(Recall)))))
	results.append([image,Trues,Detected,TP,FP,FN,Precision,Recall,Accuracy,F1])
	out.write("%s,%d,%d,%d,%d,%d,%s,%s,%s,%s\n"%(image,Trues,Detected,TP,FP,FN,Precision,Recall,Accuracy,F1))
results = np.array(results)
summ = open(sys.argv[2].split(".csv")[0]+"_accuracy_summ.csv","w")
Prec = '{:.6f}'.format(results[:,6].astype(float).sum(axis=0)/len(results[:,6]))
Reca = '{:.6f}'.format(results[:,7].astype(float).sum(axis=0)/len(results[:,7]))
Accu = '{:.6f}'.format(results[:,8].astype(float).sum(axis=0)/len(results[:,8]))
Fone = '{:.6f}'.format(results[:,9].astype(float).sum(axis=0)/len(results[:,9]))
summ.write("Precision,Recall,Accuracy,F1\n")
summ.write("%s,%s,%s,%s"%(Prec,Reca,Accu,Fone))
