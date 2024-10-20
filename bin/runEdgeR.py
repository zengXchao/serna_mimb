import argparse 
import tempfile
import os

parser = argparse.ArgumentParser(description='Generating count matrix and excuting edgeR')
parser.add_argument('-f', '--fdir', help="Directory of featureCounts results")
parser.add_argument('-o', '--outputdir', help="Directory for output files")
parser.add_argument('-t', '--treatment', help="Comma separated list of treated samples")
parser.add_argument('-c', '--control', help="Comma separated list of controlled samples)")
parser.add_argument('-g', '--gtf', help="GTF file")
args = parser.parse_args() 
controlList, treatList = [args.control], [args.treatment]
if "," in controlList[0]: controlList = controlList[0].split(",")
if "," in treatList[0]: treatList = treatList[0].split(",")

outName = "_".join(controlList[0].split("_")[:-1])+"_vs_"+"_".join(treatList[0].split("_")[:-1])
EDGER_DIR = args.outputdir

os.system("mkdir -p " + EDGER_DIR)

txList = []
with open(args.gtf) as fh:
	for line in fh:
		a = line.rstrip().split("\t")
		b = a[-1].split()
		att = {}
		for i in range(0,len(b),2):
			att[b[i]] = b[i+1][1:-2]
		if a[2] == "transcript":
			txList.append(att["transcript_id"])

counts = {}
for sample in controlList+treatList:
	counts[sample] = {}
	for tx in txList:
		counts[sample][tx] = 0
	with open(args.fdir+"/"+sample+".txt") as fh:
		fh.readline()
		fh.readline()
		for line in fh:
			a = line.rstrip().split()
			tx, val = a[0], int(a[-1])
			counts[sample][tx] = val

countFile = EDGER_DIR+"/"+outName+".count"
edgerFile = EDGER_DIR+"/"+outName+".edger"
with open(countFile,"w") as fh:
	fh.write("\t".join(["transcript_id"]+controlList+treatList)+"\n")
	for tx in txList:
		tmp = [tx]
		for sample in controlList+treatList:
			tmp.append(str(counts[sample][tx]))
		fh.write("\t".join(tmp)+"\n")

group = []
for x in controlList:
	group.append("C")
for x in treatList:
	group.append("T")

cmd = "Rscript bin/edger.R %s %s %s"%(countFile,",".join(group),edgerFile)
os.system(cmd)

countHeader = None
countData = {}
with open(countFile) as fh:
	countHeader = fh.readline().rstrip()
	for line in fh:
		a = line.rstrip().split()
		countData[a[0]] = "\t".join(a)

with open(EDGER_DIR+"/diff_exp_tx.tsv", "w") as fhW:
	fhW.write(countHeader+"\tlog2_fc\tlog10_cpm\tp_value\tfdr\n")
	with open(edgerFile) as fh:
		fh.readline()
		for line in fh:
			a = line.rstrip().split()
			fhW.write(countData[a[0].strip("\"")] + "\t" + "\t".join(a[1:]) + "\n")

