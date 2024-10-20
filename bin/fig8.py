import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib_venn import venn2
import sys

targetRNA_file = sys.argv[1]
hubRNA_file = sys.argv[2]
targetRNA_label = sys.argv[3]
hubRNA_label = sys.argv[4]
allRNA_file = sys.argv[5]
output_png = sys.argv[6]
output_pval = sys.argv[7]

def read_genes(file_path):
    with open(file_path, 'r') as file:
        genes = set(line.strip() for line in file)
    return genes

def plot_venn(set1, set2, targetRNA_label, hubRNA_label, output_file):
    plt.figure(figsize=(8, 6))
    venn= venn2([set1, set2], (targetRNA_label, hubRNA_label))
    # plt.title('Venn Diagram of Gene Lists')
    plt.savefig(output_file, dpi=100)
    plt.show()

targetRNA = read_genes(targetRNA_file)
hubRNA = read_genes(hubRNA_file)

plot_venn(targetRNA, hubRNA, targetRNA_label, hubRNA_label, output_png)

def get_gene_set(fpath):
    gene_set = set([])
    with open(fpath) as fh:
        for line in fh:
            a = line.split()
            gene_set.add(a[0].replace("_RI",""))
    return gene_set

BG = get_gene_set(allRNA_file)

n = len(BG & hubRNA)
M = len(BG) 

N = len(targetRNA)
k = len(targetRNA & hubRNA)
p = stats.hypergeom(M=M, n=n, N=N).sf(k-1)
with open(output_pval,"w") as fh:
    fh.write("%s & hubRNA: %d  (Hypergeometric test p = %s)"%(targetRNA_label, len(targetRNA & hubRNA),"{:.2e}".format(p))+"\n")
