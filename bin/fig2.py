import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import ranksums
import os
import math


seRNA_file = sys.argv[1]
exRNA_file = sys.argv[2]
apex_fpkm_dir = sys.argv[3]
output_png = sys.argv[4]
output_pval = sys.argv[5]

def getTid(fpath):
	dat = []
	with open(fpath) as fh:
		for line in fh:
			tid = line.split()[0]
			dat.append(tid)
	return dat

def getApexFpkm(fpath):
	dat = {}
	with open(fpath) as fh:
		for line in fh:
			if line[0] == "#": continue
			a = line.split()
			if a[2] != "transcript": continue
			tid, fpkm = a[-7][1:-2], float(a[-3][1:-2])
			dat[tid] = fpkm
	return dat

datalist = {"NU":{"TrackName":"Nucleus",},
	"NC":{"TrackName":"Nucleolus",},
	"NL":{"TrackName":"Nuclear lumina",},
	"NP":{"TrackName":"Nuclear pore",},
	"CY":{"TrackName":"Cytosol",},
	"EM":{"TrackName":"ER membrane cytosol facing",},
	"OM":{"TrackName":"Outer mitochondrial membrane",},
	"MM":{"TrackName":"Mitochondrial matrix",},
	"EL":{"TrackName":"ER lumen",}}

seList = getTid(seRNA_file)
exList = getTid(exRNA_file)

X, Y, Z = [],[],[]
SE,EX = {},{}
fhW = open(output_pval,"w")
for x in datalist:
	SE[x], EX[x] = [],[]
	datalist[x]["c1"] = getApexFpkm(apex_fpkm_dir+"/"+x+"_control_rep1.gtf")
	datalist[x]["c2"] = getApexFpkm(apex_fpkm_dir+"/"+x+"_control_rep2.gtf")
	datalist[x]["t1"] = getApexFpkm(apex_fpkm_dir+"/"+x+"_target_rep1.gtf")
	datalist[x]["t2"] = getApexFpkm(apex_fpkm_dir+"/"+x+"_target_rep2.gtf")

	for k in seList:
		c = (datalist[x]["c1"][k]+datalist[x]["c2"][k])/2
		t = (datalist[x]["t1"][k]+datalist[x]["t2"][k])/2
		if c < 1.0 and t < 1.0: continue 
		v = math.log((t+1)/(c+1),2)
		X.append("SE")
		Y.append(v)
		Z.append(x)
		SE[x].append(v)
	for k in exList:
		c = (datalist[x]["c1"][k]+datalist[x]["c2"][k])/2
		t = (datalist[x]["t1"][k]+datalist[x]["t2"][k])/2
		if c < 1.0 and t < 1.0: continue 
		v = math.log((t+1)/(c+1),2)
		X.append("EX")
		Y.append(v)
		Z.append(x)
		EX[x].append(v)

	w, p = ranksums(SE[x], EX[x])
	fhW.write(x + "\t" + "{:.2e}".format(p) + "\n")
fhW.close()

fig = plt.figure(figsize=(5,2.5))
ax = sns.boxplot(x=Z, y=Y, hue=X, showfliers=False, palette=["#D97A32", "#559D3E"])
ax.axhline(0, linewidth=1, color='gray', ls="--")
plt.title("Subcellular localization")
plt.ylim([-2, 3.1])
fig.tight_layout() 
fig.savefig(output_png, dpi=300)

