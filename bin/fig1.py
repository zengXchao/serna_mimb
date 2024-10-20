# import sys
# import matplotlib.pyplot as plt
# import seaborn as sns

# diff_exp_tx = sys.argv[1]
# cmp_cutoff = sys.argv[2] 
# fdr_cutoff = sys.argv[3]
# tx_fa = sys.argv[4]
# output_png = sys.argv[5]
# output_serna = sys.argv[6]
# output_exrna = sys.argv[7]

# diff_exp_tx_list = []
# with open(diff_exp_tx,"r") as fh:
# 	fh.readline()
# 	for line in fh:
# 		a = line.rstrip().split()
# 		tid, log2fc, log10cpm, fdr = a[0], float(a[-4]), float(a[-3]), float(a[-1])
# 		diff_exp_tx_list.append([tid, log2fc, log10cpm, fdr])

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gffutils
from Bio import SeqIO

diff_exp_tx = sys.argv[1]
cpm_cutoff = float(sys.argv[2])
fdr_cutoff = float(sys.argv[3])
tx_fa = sys.argv[4]
tx_gtf = sys.argv[5]
output_png = sys.argv[6]
output_serna = sys.argv[7]
output_exrna = sys.argv[8]

# 读取数据
df = pd.read_csv(diff_exp_tx, sep='\t')

# 计算log10_cpm的阈值并筛选数据
log10_cpm_cutoff = np.log10(cpm_cutoff)
df = df[df['log10_cpm'] > log10_cpm_cutoff]
df['1-fdr'] = 1 - df['fdr']

# 标记数据点颜色
df['color'] = 'grey'  # 默认颜色
up_df = df[(df['fdr'] < fdr_cutoff) & (df['log2_fc'] > 0)]
down_df = df[(df['fdr'] < fdr_cutoff) & (df['log2_fc'] < 0)]
df.loc[up_df.index, 'color'] = 'red'
df.loc[down_df.index, 'color'] = 'orange'

# 计算Non-significant数据集
remaining_indices = df[(df['color'] == 'grey')].index  # 首先获取灰色点的索引
# ns_df = df.loc[remaining_indices].nlargest(len(up_df), 'fdr')  # 按照fdr递减顺序选择
df['abs_log2_fc'] = df['log2_fc'].abs()  # 计算log2_fc的绝对值
ns_df = df.loc[remaining_indices].nsmallest(len(up_df), 'abs_log2_fc')  # 按照log2_fc绝对值递增顺序选择
df.loc[ns_df.index, 'color'] = 'blue'

# 绘图
fig, ax = plt.subplots()
colors = {'red': 'Up', 'orange': 'Down', 'blue': 'NS', 'grey': 'Others'}
# grouped = df.groupby('color')
# for key, group in grouped:
#	if key in colors:
color_order = ['grey', 'orange', 'blue', 'red']
for key in color_order:
	group = df[df['color'] == key]
	ax.scatter(group['log2_fc'], group['1-fdr'], c=key, label=f'{colors[key]} (n={len(group)})', alpha=0.6)
ax.set_xlabel('log2 fold change')
ax.set_ylabel('1 - fdr')
ax.set_title('Impr vs Conv')
ax.legend()
plt.savefig(output_png, dpi=100)

# 从FASTA文件读取序列
fasta_seqs = SeqIO.to_dict(SeqIO.parse(tx_fa, 'fasta'))

# 从GTF文件读取位置信息
db = gffutils.create_db(tx_gtf, dbfn=':memory:', force=True, merge_strategy='merge', disable_infer_transcripts=True)

# 输出特定数据
def output_data(df, output_file, fasta_seqs, db):
	with open(output_file, 'w') as f:
		for idx, row in df.iterrows():
			transcript = db[row['transcript_id']]
			sequence = fasta_seqs[row['transcript_id']].seq
			f.write(f"{row['transcript_id']}\t{transcript.chrom}\t{transcript.start}\t{transcript.end}\t{transcript.strand}\t{sequence}\n")

# 筛选符合RI条件的转录本
ri_set = {x.replace('_RI', '') for x in up_df['transcript_id'] if  x.endswith('_RI')}
base_set = {x for x in up_df['transcript_id'] if not x.endswith('_RI')}
redundant_set = ri_set.intersection(base_set)
up_df = up_df[~up_df['transcript_id'].isin(redundant_set)]

ri_set = {x.replace('_RI', '') for x in ns_df['transcript_id'] if  x.endswith('_RI')}
base_set = {x for x in ns_df['transcript_id'] if not x.endswith('_RI')}
redundant_set = ri_set.intersection(base_set)
ns_df = ns_df[~ns_df['transcript_id'].isin(redundant_set)]

# 输出Up和Non-significant数据集
print("")
print("seRNA", len(up_df))
print("exRNA", len(ns_df))
output_data(up_df, output_serna, fasta_seqs, db)
output_data(ns_df, output_exrna, fasta_seqs, db)


