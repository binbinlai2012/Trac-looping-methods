#!/usr/bin/env python

import os,sys
from string import *
import scipy
import scipy.stats

f = open(sys.argv[1], 'r')
petc = int(sys.argv[2])
thr = float(sys.argv[3])
outf = open(sys.argv[4], 'w')


#plist = []
#linelist = []
outf.write(f.readline())
#line = f.readline()
#tline = line.strip()
#outf.write(tline+"\tFDR\n")
for line in f:
	tline = line.strip()
	sline = tline.split()
#	l = int(sline[5])
#	if ( l == 2000 ):
#		continue
	pe = int(sline[9])
	
	fdr = float(sline[19])
#	if ( p <= thr ):
#	plist.append(p)
#	linelist.append(tline)
	if ( fdr <= thr and pe >= petc ):
		outf.write(line)
f.close()
outf.close()

