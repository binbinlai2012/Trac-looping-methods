#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
FDR_thr = float(sys.argv[2])
outf_inc = open(sys.argv[3], 'w')
outf_dec = open(sys.argv[4], 'w')
outf_noc = open(sys.argv[5], 'w')

header = f.readline()
outf_inc.write(header)
outf_dec.write(header)
outf_noc.write(header)
for line in f:
	tline = line.strip()
	sline = tline.split()
	fdr_BvsA = float(sline[15])
	fdr_AvsB = float(sline[16])
	norm_Fc_BvsA = float(sline[11])
	norm_Fc_AvsB = float(sline[13])
	if fdr_BvsA < FDR_thr and norm_Fc_BvsA >= 2:
		outf_inc.write(line)
	elif fdr_AvsB < FDR_thr and norm_Fc_AvsB >= 2:
		outf_dec.write(line)
	else:
		outf_noc.write(line)
f.close()
outf_inc.close()
outf_dec.close()
outf_noc.close()
