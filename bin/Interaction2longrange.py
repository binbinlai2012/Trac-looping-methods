#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
outf = open(sys.argv[2], 'w')
f.readline()

id = 0
for line in f:
	tline = line.strip()
	sline = tline.split()
	id += 1
	chr = sline[0]
	start1 = sline[1]
	end1 = sline[2]
	start2 = sline[3]
	end2 = sline[4]
	c = sline[5]
	outf.write("%s\t%s\t%s\t%s:%s-%s,%s\t%d\t.\n" % (chr, start1, end1, chr, start2, end2, c, id ))
	id+=1
	outf.write("%s\t%s\t%s\t%s:%s-%s,%s\t%d\t.\n" % (chr, start2, end2, chr, start1, end1, c, id ))
outf.close()
f.close()
