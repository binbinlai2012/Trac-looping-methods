#!/usr/bin/env python

import os,sys
from string import *
import scipy
import scipy.stats

f = open(sys.argv[1], 'r')
#thr = float(sys.argv[2])
outf = open(sys.argv[2], 'w')


plist_BvsA = []
plist_AvsB = []
linelist = []
#outf.write(f.readline())
line = f.readline()
headline = line.strip()

for line in f:
	tline = line.strip()
	sline = tline.split()
#	l = int(sline[5])
#	if ( l == 2000 ):
#		continue
	p_BvsA = float(sline[12])
	p_AvsB = float(sline[14])
#	if ( p <= thr ):
	plist_BvsA.append(p_BvsA)
	plist_AvsB.append(p_AvsB)
	linelist.append(tline)
	
f.close()
parray_BvsA = scipy.array(plist_BvsA)
prankarray_BvsA = scipy.stats.rankdata(parray_BvsA)
totalnumber_BvsA = len(plist_BvsA)

FDR_BvsA = []
for i in range(totalnumber_BvsA):
	alpha = plist_BvsA[i] * totalnumber_BvsA / prankarray_BvsA[i]
	if alpha > 1:
		alpha = 1
	FDR_BvsA.append(alpha)


parray_AvsB = scipy.array(plist_AvsB)
prankarray_AvsB = scipy.stats.rankdata(parray_AvsB)
totalnumber_AvsB = len(plist_AvsB)

FDR_AvsB = []
for i in range(totalnumber_AvsB):
	alpha = plist_AvsB[i] * totalnumber_AvsB / prankarray_AvsB[i]
	if alpha > 1:
		alpha = 1
	FDR_AvsB.append(alpha)

outf.write(headline+"\tFDR_BvsA\tFDR_AvsB\n")

for i in range(len(linelist)):
	outf.write(linelist[i]+"\t%e\t%e\n" % (FDR_BvsA[i], FDR_AvsB[i]) )
outf.close()
