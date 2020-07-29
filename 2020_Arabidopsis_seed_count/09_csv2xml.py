#! /usr/bin/env python

import csv
import sys,os
inp = sys.argv[1]
csv.register_dialect('custom',
                     delimiter=',',
                     doublequote=True,
                     escapechar=None,
                     quotechar='"',
                     quoting=csv.QUOTE_MINIMAL,
                     skipinitialspace=False)

oneline = open(inp,"r").readline()
width = oneline.strip().split(",")[1]
height = oneline.strip().split(",")[2]
with open(inp) as ifile:
	data = csv.reader(ifile, dialect='custom')
	name = inp.split(".csv")[0]
	print "<annotation>"
	print "<folder>annotation</folder>"
	print "<filename>"+name+"</filename>"
	print "<path>annotaion/"+name+"</path>"
	print "<source>"
    	print "<database>Unknown</database>"
  	print "</source>"
	print "<size>"
    	print "<width>"+width+"</width>"
    	print "<height>"+height+"</height>"
    	print "<depth>3</depth>"
  	print "</size>"
  	print "<segmented>0</segmented>"
    	for record in data:
		print "<object>"
    		print "<name>seed</name>"
    		print "<pose>Unspecified</pose>"
    		print "<truncated>0</truncated>"
    		print "<difficult>0</difficult>"
    		print "<bndbox>"
      		print "<xmin>"+record[4]+"</xmin>"
      		print "<ymin>"+record[5]+"</ymin>"
      		print "<xmax>"+record[6]+"</xmax>"
      		print "<ymax>"+record[7]+"</ymax>"
    		print "</bndbox>"
  		print "</object>"
    	print "</annotation>"
