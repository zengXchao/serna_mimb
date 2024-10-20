import sys
import statistics
import math

stringtie_dir = sys.argv[1]
sample_list = sys.argv[2:]

meanThreshold = 0.1

fpkm = {}
for x in sample_list:
	with open("%s/%s.gtf"%(stringtie_dir,x)) as fh:
		for line in fh:
			if line[0] == "#": continue
			a = line.split()
			if a[2] != "transcript": continue
			fid, val = a[11][1:-2], float(a[15][1:-2])
			if fid not in fpkm: fpkm[fid] = []
			fpkm[fid].append(val)

## FPKM of exons and introns
fpkmMeanExon, fpkmMeanIntron = {},{}
fpkmMeanExonFilteredLog, fpkmMeanIntronFilteredLog = {},{}
fh = open("%s/fpkm_of_exon_intron.txt"%stringtie_dir,"w") 
for fid in fpkm:
	m = statistics.mean(fpkm[fid])
	if "_exon_" in fid:
		fpkmMeanExon[fid] = m
		if m >= meanThreshold: 
			fpkmMeanExonFilteredLog[fid] = math.log(m,10)
	elif "_intron_" in fid:
		fpkmMeanIntron[fid] = m
		if m >= meanThreshold: 
			fpkmMeanIntronFilteredLog[fid] = math.log(m,10)
	fh.write(fid+"\t"+str(m)+"\n")
fh.close()

## Retaintion score (RS) of introns
rsIntron = {}
fh = open("%s/intron_rs.txt"%stringtie_dir,"w")

fh.write("# Intron FPKM\n")
fh.write("# "+"\t".join(["intronID","intronFPKM","flankingExonFPKM","RS"])+"\n")
for fid in fpkmMeanIntron:
	a = fid.split("_intron_")
	intronFPKM = fpkmMeanIntron[fid]
	preExonFPKM = 0
	if a[0]+"_exon_"+a[1] in fpkmMeanExon:
		fpkmMeanExon[a[0]+"_exon_"+a[1]]
	nextExonFPKM = 0
	if a[0]+"_exon_"+str(int(a[1])+1) in fpkmMeanExon:
		nextExonFPKM = fpkmMeanExon[a[0]+"_exon_"+str(int(a[1])+1)]
	exonFPKM = (preExonFPKM+nextExonFPKM)/2.
	rs = 0
	if exonFPKM>0: rs = intronFPKM/exonFPKM
	fh.write("\t".join([fid, str(intronFPKM),str(exonFPKM),str(rs)])+"\n")
	if intronFPKM > 0 and intronFPKM >= 1 and exonFPKM > 0 and exonFPKM >= 1: 
		rsIntron[fid] = rs
fh.close()