import os
from scipy import stats
import math
import sys

seRNA_file = sys.argv[1]
exRNA_file = sys.argv[2]
meme_db_file = sys.argv[3]
motif_tmp_dir = sys.argv[4]

os.system("mkdir -p " + motif_tmp_dir+"/fimo/seq/se")
os.system("mkdir -p " + motif_tmp_dir+"/fimo/seq/ex")

with open(motif_tmp_dir+"/fimo/SE.fa", "w") as fh:
	with open(seRNA_file) as fh2:
		for line in fh2:
			[tid,seq] = line.rstrip().split()
			fh.write(">"+tid+"\n"+seq+"\n")
			seqFile = motif_tmp_dir+"/fimo/seq/se/"+tid+".fa"
			fimoDir = motif_tmp_dir+"/fimo/seq/se/"+tid
			with open(seqFile, "w") as fh2:
				fh2.write(">"+tid+"\n"+seq+"\n")
			os.system("fimo --verbosity 1 --thresh 0.01 --norc --o %s %s %s"%(fimoDir, meme_db_file, seqFile))

with open(motif_tmp_dir+"/fimo/EX.fa", "w") as fh:
	with open(exRNA_file) as fh2:
		for line in fh2:
			[tid,seq] = line.rstrip().split()
			fh.write(">"+tid+"\n"+seq+"\n")
			seqFile = motif_tmp_dir+"/fimo/seq/ex/"+tid+".fa"
			fimoDir = motif_tmp_dir+"/fimo/seq/ex/"+tid
			with open(seqFile, "w") as fh2:
				fh2.write(">"+tid+"\n"+seq+"\n")
			os.system("fimo --verbosity 1 --thresh 0.01 --norc --o %s %s %s"%(fimoDir, meme_db_file, seqFile))

os.system("cat " + motif_tmp_dir+"/fimo/seq/se/*/fimo.tsv > " + motif_tmp_dir+"/fimo/SE.fimo")
os.system("cat " + motif_tmp_dir+"/fimo/seq/ex/*/fimo.tsv > " + motif_tmp_dir+"/fimo/EX.fimo")

## Removed
## os.system("fimo --verbosity 1 --norc --thresh 0.01 --qv-thresh --max-stored-scores 1000000000000 --o %s %s %s"%(motif_tmp_dir+"/fimo/se", motif_tmp_dir+"/fimo/CISBP-RNA_v0.6_RBP.motif.meme", motif_tmp_dir+"/fimo/SE.fa"))
## os.system("fimo --verbosity 1 --norc --thresh 0.01 --qv-thresh --max-stored-scores 1000000000000 --o %s %s %s"%(motif_tmp_dir+"/fimo/ex", motif_tmp_dir+"/fimo/CISBP-RNA_v0.6_RBP.motif.meme", motif_tmp_dir+"/fimo/EX.fa"))	

binN = 5
def load_tx_len(fpath):
	txLen = {}
	with open(fpath) as fh:
		for line in fh:
			[tid,seq] = line.rstrip().split()
			val = len(seq)
			if val < 100: continue
			txLen[tid] = val
	return txLen

seLen = load_tx_len(seRNA_file)
exLen = load_tx_len(exRNA_file)

def load_fimo_data(fpath, sLen):
	df1, df2 = {},{}
	dat1,dat2 = {},{}
	with open(fpath) as fh:
		for line in fh:
			if line[0] != "M": continue
			a = line.split()
			sId = a[2]
			if sId not in sLen: continue
			motif, start = a[1], int(a[3])
			if motif not in df1:
				df1[motif] = {}
				df2[motif] = {}
				for x in sLen:
					df1[motif][x] = [0.] * binN
					df2[motif][x] = 0.
			binIdx = math.floor((start-1) / math.ceil(sLen[sId]/binN))
			df1[motif][sId][binIdx] += 1
			df2[motif][sId] += 1
	for m in df2:
		dat2[m] = 0.
		for k in df2[m]:
			dat2[m] += df2[m][k]/sLen[k]
		dat2[m] /= len(sLen)
	for m in df1:
		dat1[m] = [0.] * binN
		for k in df1[m]:
			for i in range(binN):
				dat1[m][i] += df1[m][k][i]/math.ceil(sLen[k]/binN)
		for i in range(binN):
			dat1[m][i] /= len(sLen)
	return dat1, dat2

seDat1, seDat2 = load_fimo_data(motif_tmp_dir+"/fimo/SE.fimo", seLen)
exDat1, exDat2 = load_fimo_data(motif_tmp_dir+"/fimo/EX.fimo", exLen)

def write_dat(fpath,dat1,dat2):
	with open(fpath,"w") as fh:
		for k in dat1:
			tmp = [k]
			tmp.extend([str(round(x,4)) for x in dat1[k]])
			tmp.append(str(round(dat2[k],4)))
			fh.write("\t".join(tmp)+"\n")
write_dat(motif_tmp_dir+"/fimo/SE.txt",seDat1, seDat2)
write_dat(motif_tmp_dir+"/fimo/EX.txt",exDat1, exDat2)

def write_dat2(fpath,dat1,dat2,dat3,dat4):
	with open(fpath,"w") as fh:
		for k in dat1:
			if k not in dat3: continue
			tmp = [k]
			for i in range(binN):
				try:
					tmp.append(str(round(math.log((dat1[k][i])/(dat3[k][i]),2),3)))
				except Exception as e:
					tmp.append("0")
			try:
				tmp.append(str(round(math.log((dat2[k])/(dat4[k]),2),3)))
			except Exception as e:
				tmp.append("0")
			
			fh.write("\t".join(tmp)+"\n")
write_dat2(motif_tmp_dir+"/RBP_preference.txt",seDat1, seDat2, exDat1, exDat2)
