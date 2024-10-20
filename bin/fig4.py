import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import ranksums
import sys
import os

merged_fa = sys.argv[1]
seRNA_file = sys.argv[2]
exRNA_file = sys.argv[3]
output_tmp = sys.argv[4]
output_png = sys.argv[5]
output_pval = sys.argv[6]

def load_tid(fpath):
	dat = {}
	with open(fpath) as fh:
		for line in fh:
			a = line.split()
			dat[a[0]] = 0
	return dat


all_seq = open(output_tmp + "/all_seq.tsv","w")
se_seq = open(output_tmp + "/se_seq.tsv","w")
ex_seq = open(output_tmp + "/ex_seq.tsv","w")
se_list = load_tid(seRNA_file)
ex_list = load_tid(exRNA_file)
seq_name, seq = None, None
with open(merged_fa) as fh:
	for line in fh:
		if line[0] == ">":
			if seq_name is not None:
				all_seq.write(seq_name+"\t"+seq+"\n")
				if seq_name in se_list:
					se_seq.write(seq_name+"\t"+seq+"\n")
				if seq_name in ex_list:
					ex_seq.write(seq_name+"\t"+seq+"\n")
			seq_name = line.rstrip()[1:]
			seq = ""
		else:
			seq += line.rstrip()
all_seq.write(seq_name+"\t"+seq+"\n")
if seq_name in se_list:
	se_seq.write(seq_name+"\t"+seq+"\n")
if seq_name in ex_list:
	ex_seq.write(seq_name+"\t"+seq+"\n")
all_seq.close()
se_seq.close()
ex_seq.close()

def generate_terminal_fasta(seqFile,prefix,l):
	fh5 = open(prefix+"_5end.fa","w")
	fh3 = open(prefix+"_3end.fa","w")
	head = None
	with open(seqFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			head = ">" + a[0]
			seq = a[1]
			seqLen = len(seq)
			if seqLen < 2*l: continue
			fh5.write(head+"\n")
			fh3.write(head+"\n")
			fh5.write(seq[:l]+"\n")
			fh3.write(seq[-l:]+"\n")
	fh5.close()
	fh3.close()

def get_mfe(fipath, fopath):
	dat = []
	with open(fopath,"w") as fhO:
		tx = None
		with open(fipath) as fhI:
			for line in fhI:
				if line[0] == ">":
					tx = line[1:].split(":")[0]
				elif line[-2] == ")":
					mfe = line[:-2].split("(")[-1]
					fhO.write(tx+"\t"+mfe+"\n")
					dat.append(float(mfe))
	return dat

def run_RNAfold():
	generate_terminal_fasta(output_tmp + "/se_seq.tsv",output_tmp + "/SE",300)
	generate_terminal_fasta(output_tmp + "/ex_seq.tsv",output_tmp + "/EX",300)
	generate_terminal_fasta(output_tmp + "/all_seq.tsv",output_tmp + "/BG",300)
	os.system("RNAfold --noPS -i %s > %s"%(output_tmp + "/SE_5end.fa", output_tmp + "/SE_5end.out"))
	os.system("RNAfold --noPS -i %s > %s"%(output_tmp + "/SE_3end.fa", output_tmp + "/SE_3end.out"))
	os.system("RNAfold --noPS -i %s > %s"%(output_tmp + "/EX_5end.fa", output_tmp + "/EX_5end.out"))
	os.system("RNAfold --noPS -i %s > %s"%(output_tmp + "/EX_3end.fa", output_tmp + "/EX_3end.out"))
	os.system("RNAfold --noPS -i %s > %s"%(output_tmp + "/BG_5end.fa", output_tmp + "/BG_5end.out"))
	os.system("RNAfold --noPS -i %s > %s"%(output_tmp + "/BG_3end.fa", output_tmp + "/BG_3end.out"))


def plotMFE():
	run_RNAfold()

	SE5 = get_mfe(output_tmp + "/SE_5end.out", output_tmp + "/SE_5end.mfe")
	SE3 = get_mfe(output_tmp + "/SE_3end.out", output_tmp + "/SE_3end.mfe")
	EX5 = get_mfe(output_tmp + "/EX_5end.out", output_tmp + "/EX_5end.mfe")
	EX3 = get_mfe(output_tmp + "/EX_3end.out", output_tmp + "/EX_3end.mfe")
	BG5 = get_mfe(output_tmp + "/BG_5end.out", output_tmp + "/BG_5end.mfe")
	BG3 = get_mfe(output_tmp + "/BG_3end.out", output_tmp + "/BG_3end.mfe")

	X,Y,Z = [],[],[]
	for v in SE5:
		X.append("5' end")
		Y.append(v)
		Z.append("SE")
	for v in EX5:
		X.append("5' end")
		Y.append(v)
		Z.append("EX")
	for v in BG5:
		X.append("5' end")
		Y.append(v)
		Z.append("BG")
	for v in SE3:
		X.append("3' end")
		Y.append(v)
		Z.append("SE")
	for v in EX3:
		X.append("3' end")
		Y.append(v)
		Z.append("EX")
	for v in BG3:
		X.append("3' end")
		Y.append(v)
		Z.append("BG")

	fig = plt.figure(figsize=(3,2.5))
	ax = sns.boxplot(x=X, y=Y, hue=Z,  showfliers=False, palette=["#D97A32", "#559D3E", "#555555"])
	plt.title("MFE")
	y1,y2 = ax.get_ylim()
	plt.ylim([y1, 1.5*300/10])
	plt.legend(loc=4)
	fig.tight_layout() 
	fig.savefig(output_png, dpi=300)

	fh = open(output_pval,"w")
	w, p = ranksums(SE5, EX5)
	fh.write("5end (300bp):  SE vs EX Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	w, p = ranksums(SE5, BG5)
	fh.write("5end (300bp):  SE vs BG Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	w, p = ranksums(SE3, EX3)
	fh.write("3end (300bp):  SE vs EX Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	w, p = ranksums(SE3, BG3)
	fh.write("3end (300bp):  SE vs BG Wilcoxon rank-sum test p = " + "{:.2e}".format(p) + "\n")
	fh.close()

plotMFE()
