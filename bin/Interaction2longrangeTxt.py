#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
outf = open(sys.argv[2], 'w')
f.readline()


for line in f:
	tline = line.strip()
	sline = tline.split()
	
	chr = sline[0]
	start1 = sline[1]
	end1 = sline[2]
	start2 = sline[3]
	end2 = sline[4]
	c = sline[5]
	outf.write("%s,%s,%s\t%s,%s,%s\t%s\n" % (chr, start1, end1, chr, start2, end2, c ))
	
outf.close()
f.close()
