#!/usr/bin/env python

import os,sys
from string import *

bedpefile = sys.argv[1]
out150file = sys.argv[2]
out1K2Mfile = sys.argv[3]
outstatfile = sys.argv[4]

interc = 0
intrac = 0
chrM = 0
intrac_150 = 0
intrac_1k = 0
intrac_200k = 0
intrac_200kl = 0

cm1k = 0
cm0 = 0

f=open(bedpefile, 'r')
outf1 = open(out150file, 'w')
outf2 = open(out1K2Mfile, 'w')
for line in f:
	tline = line.strip()
	sline = tline.split()
	chr1 = sline[0]
	chr2 = sline[3]
	
	if chr1 == "chrM" or chr2 == "chrM":
		chrM += 1
		continue
	elif chr1 != chr2:
		interc += 1
	
	
	start1 = int(sline[1])
	end1 = int(sline[2])
	
	start2 = int(sline[4])
	end2 = int(sline[5])
	str1 = sline[7]
	str2 = sline[8]
	sc = int(sline[6])
	c1 = start1
	if str == '-':
		c1 = end1
	c2 = start2
	if str == '-':
		c2 = end2
		
	l = c2 - c1
	if l <= 0:
		if l < -1000:
			cm1k += 1
		else:
			cm0 += 1
		continue
	elif l < 150:
		outf1.write("%s\t%d\t%d\t%d\t%d\t+\n" % (chr1, c1, c2, l, sc ))
		intrac_150 += 1
	elif l < 1000:
		intrac_1k += 1
	else:
		if l < 2000000:
			outf2.write(line)
		
		if l < 200000:
			intrac_200k += 1
		else:
			intrac_200kl += 1
f.close()
outf1.close()
outf2.close()

intrac = intrac_150 + intrac_1k + intrac_200k + intrac_200kl

outfs = open(outstatfile, 'w')
outfs.write("#chrM_PETs\t%d\n" % chrM)
outfs.write("#interChr_PETs\t%d\n" % interc)
outfs.write("#intraChr_PETs\t%d\n" % intrac)
outfs.write("#intraChr_PETs(<150bp)\t%d\n" % intrac_150)
outfs.write("#intraChr_PETs(150bp,1kbp)\t%d\n" % intrac_1k)
outfs.write("#intraChr_PETs(1kbp,200kbp)\t%d\n" % intrac_200k)
outfs.write("#intraChr_PETs(>200kbp)\t%d\n" % intrac_200kl)
outfs.write("##intraChr_PETs(-1k,0bp)\t%d\n" % cm0)
outfs.write("##intraChr_PETs(<-1kbp)\t%d\n" % cm1k)


			
		
		
		
		
		
		
		
		
