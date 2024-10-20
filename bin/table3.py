import networkx as nx
import pandas as pd
import statistics
import community as community_louvain
import sys

input_csv = sys.argv[1]
output_cluster = sys.argv[2]
edge_weight_threshold = float(sys.argv[3])
seRNA_file = sys.argv[4]
exRNA_file = sys.argv[5]
motif_tmp_dir = sys.argv[6]
output_rbp = sys.argv[7]

# # Read the CSV file
# df = pd.read_csv(input_csv, index_col=0)

# min_value = df.min().min()
# max_value = df.max().max()
# df_normalized = (df - min_value) / (max_value - min_value)


# # Create a graph
# G = nx.Graph()

# # Add nodes
# for node in df.index:
#     G.add_node(node)

# # Add edges with weights
# for i in df_normalized.index:
#     for j in df_normalized.columns:
#         if i != j and not pd.isnull(df_normalized.loc[i, j]) and df_normalized.loc[i, j] >= edge_weight_threshold:
#             G.add_edge(i, j, weight=df_normalized.loc[i, j])

# # Run the Louvain algorithm for community detection
# partition = community_louvain.best_partition(G, weight='weight')

# # Print nodes and their community labels
# with open(output_cluster,"w") as fh:
#     for node, community in partition.items():
#         community += 1
#         fh.write(f'{node}\t{community}\n')

clusterList={}
with open(output_cluster) as fh:
    fh.readline()
    for line in fh:
        a = line[1:-1].split()
        if a[1] not in clusterList:
            clusterList[a[1]] = []
        clusterList[a[1]].append(a[0])

def load_tx_len(fpath, whitelist=None):
    txLen = {}
    with open(fpath) as fh:
        for line in fh:
            [tid,seq] = line.rstrip().split()
            val = len(seq)
            if val < 100: continue
            if whitelist != None and tid not in whitelist: continue
            txLen[tid] = val
    return txLen

seLen = {}
for c in clusterList:
    seLen[c] = load_tx_len(seRNA_file, clusterList[c])
nsLen = load_tx_len(exRNA_file)

def load_fimo_data(fpath, txLen):
    rbpDense = {}
    with open(fpath) as fh:
        for line in fh:
            if line[0] != "M": continue
            a = line.split()
            tx = a[2].split("::")[0]
            rbp = a[1]
            if tx not in txLen: continue
            if rbp not in rbpDense:
                rbpDense[rbp] = {}
                for tx in txLen:
                    rbpDense[rbp][tx] = 0
            rbpDense[rbp][tx] += 1
    dat = {}
    for rbp in rbpDense:
        dat[rbp] = []
        for tx in txLen:
            rbpDense[rbp][tx] = float(rbpDense[rbp][tx])/txLen[tx]
            dat[rbp].append(rbpDense[rbp][tx])
    return dat

rbpDense = {}
for c in clusterList:
    rbpDense[c] = load_fimo_data(motif_tmp_dir+"/fimo/SE.fimo",seLen[c])
rbpDense["EX"] = load_fimo_data(motif_tmp_dir+"/fimo/EX.fimo",nsLen)


rbpList = set([])
for c in rbpDense:
    for x in rbpDense[c]:
        rbpList.add(x)

# with open(FIG_4_DIR+"/k%s_d%s_rbp_profile_mean.tsv"%(kmer,d),"w") as fh:
#     fh.write("\t".join([""] + sorted(clusterList)+["EX"])+"\n")
#     for rbp in sorted(rbpList):
#         tmp = [rbp]
#         for c in sorted(clusterList)+["EX"]:
#             try:
#                 tmp.append(str(statistics.mean(rbpDense[c][rbp])))
#             except Exception as e:
#                 tmp.append("0")
#         fh.write("\t".join(tmp)+"\n")


with open(output_rbp,"w") as fh:
    fh.write("\t".join(["RBP"] + [ "cluster_"+x for x in sorted(clusterList)]+["EX"])+"\n")
    for rbp in sorted(rbpList):
        tmp = [rbp]
        for c in sorted(clusterList)+["EX"]:
            try:
                tmp.append(str(statistics.median(rbpDense[c][rbp])))
            except Exception as e:
                tmp.append("0")
        fh.write("\t".join(tmp)+"\n")