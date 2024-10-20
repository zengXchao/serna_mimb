import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import ranksums
import sys

tx_fai = sys.argv[1]
seRNA_file = sys.argv[2]
exRNA_file = sys.argv[3]
output_png = sys.argv[4]
output_pval = sys.argv[5]

txLen = {}
with open(tx_fai) as fh:
	for line in fh:
		a = line.split()
		txLen[a[0]] = float(a[1])

def load_tid(fpath):
	dat = {}
	with open(fpath) as fh:
		for line in fh:
			a = line.split()
			dat[a[0]] = 0
	return dat

seList = load_tid(seRNA_file)
exList = load_tid(exRNA_file)


SE = [txLen[x] for x in seList]
EX = [txLen[x] for x in exList]
BG = [txLen[x] for x in txLen]
X,Y,Z = [],[],[]
for val in SE:
	X.append("SE")
	Y.append(val)
	Z.append(["SE",val])
for val in EX:
	X.append("EX")
	Y.append(val)
	Z.append(["EX",val])
for val in BG:
	X.append("BG")
	Y.append(val)
	Z.append(["BG",val])

df = pd.DataFrame(Z,columns=["Category","Length"])
fig = plt.figure(figsize=(3,3))
sns.ecdfplot(df,x="Length",hue="Category",palette=["#D97A32","#559D3E","#555555"])
plt.xscale("log")
fig.tight_layout() 
fig.savefig(output_png, dpi=300)

fh = open(output_pval,"w")
w, p = ranksums(SE, EX)
fh.write("length (SE vs EX): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
w, p = ranksums(SE, BG)
fh.write("length (SE vs BG): Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
fh.close()