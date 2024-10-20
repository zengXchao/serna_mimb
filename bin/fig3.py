import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import ranksums
import sys

tx_fai = sys.argv[1]
se_bed = sys.argv[2]
ex_bed = sys.argv[3]
bg_bed = sys.argv[4]
se_x_imargi = sys.argv[5]
ex_x_imargi = sys.argv[6]
bg_x_imargi = sys.argv[7]
output_png = sys.argv[8]
output_pval = sys.argv[9]

def get_chromatin_interaction_density(bedFile, xFile, dname, txLen, Z):
	V,Dense = [],{}
	with open(bedFile) as fh:
		for line in fh:
			a = line.split()
			tid = a[3].split(":")[0]
			Dense[tid] = 0
	with open(xFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			tid = a[3].split(":")[0]
			bpN = float(a[-1])
			Dense[tid] += bpN/txLen[tid]
	for k in Dense:
		V.append(Dense[k])
		Z.append([dname,Dense[k]])
	return V 

txLen  = {}
with open(tx_fai) as fh:
	for line in fh:
		a = line.split()
		txLen[a[0]] = float(a[1])

Z = []
SE = get_chromatin_interaction_density(se_bed, se_x_imargi, "SE", txLen, Z)
EX = get_chromatin_interaction_density(ex_bed, ex_x_imargi, "EX", txLen, Z)
BG = get_chromatin_interaction_density(bg_bed, bg_x_imargi, "BG", txLen, Z)

df = pd.DataFrame(Z,columns=["Category","Chromatin interaction per nt"])
png = output_png
fig = plt.figure(figsize=(3,3))
sns.ecdfplot(df,x="Chromatin interaction per nt",hue="Category",palette=["#D97A32","#559D3E","#555555"])
plt.xlim([0, 20])
fig.tight_layout() 
fig.savefig(png, dpi=300)

fh = open(output_pval,"w")
w, p = ranksums(SE, EX)
fh.write("Chromatin interaction (SE vs EX): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
w, p = ranksums(SE, BG)
fh.write("Chromatin interaction (SE vs BG): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
w, p = ranksums(EX, BG)
fh.write("Chromatin interaction (EX vs BG): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
fh.close()
