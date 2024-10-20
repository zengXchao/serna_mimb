import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
from scipy.stats import ranksums
import sys
import os

txSeq_file = sys.argv[1]
seSeq_file = sys.argv[2]
exSeq_file = sys.argv[3]
output_png = sys.argv[4]
output_pval = sys.argv[5]

txSeq, txLen, seList, exList = {}, {}, [], []
with open(txSeq_file) as fh:
	for line in fh:
		a = line.rstrip().split()
		txSeq[a[0]] = a[1]
		txLen[a[0]] = len(a[1])
with open(seSeq_file) as fh:
	for line in fh:
		a = line.rstrip().split()
		seList.append(a[0])
with open(exSeq_file) as fh:
	for line in fh:
		a = line.rstrip().split()
		exList.append(a[0])

SE, EX, BG = [], [], []
X,Y,Z = [],[],[]
for tid in seList:
	val = 100*(txSeq[tid].count("G")+txSeq[tid].count("C"))/float(txLen[tid])
	X.append("SE")
	Y.append(val)
	Z.append(["SE",val])
	SE.append(val)
for tid in exList:
	val = 100*(txSeq[tid].count("G")+txSeq[tid].count("C"))/float(txLen[tid])
	X.append("EX")
	Y.append(val)
	Z.append(["EX",val])
	EX.append(val)
for tid in txSeq:
	val = 100*(txSeq[tid].count("G")+txSeq[tid].count("C"))/float(txLen[tid])
	X.append("BG")
	Y.append(val)
	Z.append(["BG",val])
	BG.append(val)


df = pd.DataFrame(Z,columns=["Category","GC-content (%)"])
fig = plt.figure(figsize=(3,3))
sns.ecdfplot(df,x="GC-content (%)",hue="Category",palette=["#D97A32","#559D3E","#555555"])
fig.tight_layout() 
fig.savefig(output_png, dpi=300)

fh = open(output_pval,"w")
w, p = ranksums(SE, EX)
fh.write("GC-content (SE vs EX): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
w, p = ranksums(SE, BG)
fh.write("GC-content (SE vs BG): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
fh.close()