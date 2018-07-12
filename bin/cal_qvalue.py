#!/usr/bin/env python

import os,sys
from string import *
import scipy
import scipy.stats

f = open(sys.argv[1], 'r')
#thr = float(sys.argv[2])
outf = open(sys.argv[2], 'w')


plist = []
linelist = []
#outf.write(f.readline())
line = f.readline()
tline = line.strip()
outf.write(tline+"\tFDR\n")
for line in f:
	tline = line.strip()
	sline = tline.split()
#	l = int(sline[5])
#	if ( l == 2000 ):
#		continue
	p = float(sline[18])
#	if ( p <= thr ):
	plist.append(p)
	linelist.append(tline)
f.close()
parray = scipy.array(plist)
prankarray = scipy.stats.rankdata(parray)
totalnumber = len(plist)

for i in range(totalnumber):
	alpha = plist[i] * totalnumber / prankarray[i]
	if alpha > 1:
		alpha = 1
#	if alpha < thr:
	outf.write(linelist[i]+"\t%e\n" % alpha)
outf.close()
