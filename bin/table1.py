#!/usr/bin/envi python
import os
import sys

SE_BED = sys.argv[1]
EX_BED = sys.argv[2]
BG_BED = sys.argv[3]
SE_x_chromHMM = sys.argv[4]
EX_x_chromHMM = sys.argv[5]
BG_x_chromHMM = sys.argv[6]
OUTPUT_FILE = sys.argv[7]

def get_bedLen(fpath):
	totalLen = 0
	with open(fpath) as fh:
		for line in fh:
			a = line.split()
			totalLen += int(a[2]) - int(a[1]) + 1
	return totalLen

seLen,exLen,bgLen=get_bedLen(SE_BED),get_bedLen(EX_BED),get_bedLen(BG_BED)

seg2cluster = {"Tss":0,"TssF":0,
	"PromF":1,
	"PromP":2,
	"Enh":3,"EnhF":3,
	"EnhWF":4,"EnhW":4,"DnaseU":4,"DnaseD":4,"FaireW":4,
	"CtcfO":5,"Ctcf":5,
	"Gen5\'":6,"Elon":6,"ElonW":6,"Gen3\'":6,"Pol2":6,"H4K20":6,
	"Low":7,
	"ReprD":8,"Repr":8,"ReprW":8,
	"Quies":9,"Art":9,}

cluster2annotation = ["Active Promoter","Promoter Flanking","Inactive Promoter","Candidate Strong enhancer","Candidate Weak enhancer/DNase",
	"Distal CTCF/Candidate Insulator","Transcription associated","Low activity proximal to active states","Polycomb repressed",
	"Heterochromatin/Repetitive/Copy Number Variation",]

def cal_seg_frac(totalLen, fpath):
	frac = [0.] * 10
	with open(fpath) as fh:
		for line in fh:
			a = line.rstrip().split("\t")
			seg, overlap = a[3], int(a[-1])
			frac[seg2cluster[seg]] += overlap
	for i in range(10):
		frac[i] = 100*frac[i]/totalLen
	return frac

seFrac = cal_seg_frac(seLen, SE_x_chromHMM)
exFrac = cal_seg_frac(exLen, EX_x_chromHMM)
bgFrac = cal_seg_frac(bgLen, BG_x_chromHMM)

def str2(val):
	return str(round(val,3))

with open(OUTPUT_FILE,"w") as fh:
	fh.write("Chromatin states\t%SE\t%EX\t%BG\tSE/EX\tSE/BG\n")
	for i in range(10):
		tmp = [cluster2annotation[i],str2(seFrac[i]),str2(exFrac[i]),str2(bgFrac[i])]
		try:
			tmp.append(str2(round(seFrac[i]/exFrac[i],2)))
		except:
			tmp.append("NA")
		try:
			tmp.append(str2(round(seFrac[i]/bgFrac[i],2)))
		except:
			tmp.append("NA")
		fh.write("\t".join(tmp)+"\n")









