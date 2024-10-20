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
se_x_te = sys.argv[5]
ex_x_te = sys.argv[6]
bg_x_te = sys.argv[7]
output_png = sys.argv[8]
output_pval = sys.argv[9]

def get_te_density(bedFile, xFile, dname, X, Y, Z, txLen):
	teF_set = set([])
	with open(xFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			b = a[9].split("/")
			teF = b[0]
			if "?" not in teF and "Unknown" not in teF: 
				if teF not in ["RC","RNA","Retroposon","Satellite","rRNA","scRNA","sbRNA","srpRNA","tRNA", "snRNA", "Simple_repeat", "Low_complexity"]:
					teF_set.add(teF)
	V,Dense = {},{}
	for k in teF_set:
		V[k] = []
		Dense[k] = {}
	with open(bedFile) as fh:
		for line in fh:
			a = line.split()
			for x in teF_set:
				Dense[x][a[3].split(":")[0]] = 0
	with open(xFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			b = a[9].split("/")
			teF = b[0]
			if teF in teF_set:
				bpN = float(a[-1])
				Dense[teF][a[3].split(":")[0]] += bpN/txLen[a[3].split(":")[0]]
	for k in teF_set: 
		for m in Dense[k]:
			X.append(dname)
			Y.append(Dense[k][m])
			Z.append(k)
			V[k].append(Dense[k][m])
	return V  


txLen = {}
with open(tx_fai) as fh:
	for line in fh:
		a = line.split()
		txLen[a[0]] = float(a[1])

X,Y,Z = [],[],[]
SE = get_te_density(se_bed, se_x_te, "SE", X, Y, Z, txLen)
NS = get_te_density(ex_bed, ex_x_te, "NS", X, Y, Z, txLen)
BG = get_te_density(bg_bed, bg_x_te, "BG", X, Y, Z, txLen)

fig = plt.figure(figsize=(5,2.5))
sns.boxplot(x=Z, y=Y, hue=X, showfliers=False, palette=["#D97A32", "#559D3E", "#555555"])
plt.title("Repeat density")
plt.ylim([-0.05, 0.7])
fig.tight_layout() 
fig.savefig(output_png, dpi=300)

fh = open(output_pval,"w")
for k in SE:
	fh.write(k + "\n")
	w, p = ranksums(SE[k], NS[k])
	fh.write("TE density (SE vs NS): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	w, p = ranksums(SE[k], BG[k])
	fh.write("TE density (SE vs BG): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	w, p = ranksums(NS[k], BG[k])
	fh.write("TE density (EX vs BG): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	fh.write("\n")
fh.close()
